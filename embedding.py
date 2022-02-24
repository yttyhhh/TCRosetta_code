import pandas as pd
import numpy as np
import scprep
import base64
import json
import os
import tempfile
import scipy.stats as stats
# from scipy.cluster.vq import kmeans2
# from sklearn.preprocessing import StandardScaler
from sklearn import preprocessing
# from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
# import umap
# import matplotlib
# from sklearn.preprocessing import MinMaxScaler
# from scipy.cluster.vq import kmeans2
from functools import reduce
from collections import Counter
# import sys
import pickle
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import normal_data_processing
# import sklearn.decomposition
#from Example import mycol

embedding_api = Blueprint('embedding', __name__)

# background data
background = pd.read_hdf("/workspace2/yuet/tcrosetta/background_data/background_motif.hdf")
#disease data
covid =pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/covid.hdf")
nsclc =pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/nsclc.hdf")
mel =pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/mel.hdf")
breast = pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/breast.hdf")
cmv = pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/cmv.hdf")
chy= pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/classical_Hodgkin_lymphoma.hdf")
cro=pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/Crohn_disease.hdf")
yel=pd.read_hdf("/workspace2/yuet/tcrosetta/phate_data/tcrosetta/08_filter_by_background/Yellow_fever_vaccine.hdf")
#phate_operator and phate_data
f_name = '/workspace/yuet/tcr_tool/Assesst/traind_disease_phate.sav'
filename = '/workspace/yuet/tcr_tool/Assesst/trained_phate.npy'
loaded_model = pickle.load((open(f_name, 'rb')))
load_phate = np.load(filename)

def splitSeqs(seq,count=1):
    test = Counter([ seq[pos:pos+3] for pos in range(1, len(seq)-2)])
    return Counter({key:value*count for key,value in test.items()})

def transformAAseq(ll):
    return reduce(lambda x,y : x+y , [splitSeqs(seq,ct) for (seq,ct) in ll ])

def pre_embedding_data(tab):
    tab = tab[["AASeq","Ratio"]]
    #tab = pd.read_csv(path,usecols=["AASeq","cloneFraction"])
    df = zip(tab["AASeq"],tab["Ratio"])
    res = transformAAseq(df)
    td = pd.DataFrame()
    td["AA"] = res.keys()
    td["Freq"] = res.values()
    td["ID"] = "InSample"
    td["Condition"] = 'New'
    s_pivot = pd.pivot_table(td,index="ID",columns="AA",values="Freq")
    s_pivot = s_pivot.fillna(0)
    return s_pivot

def filter_by_librarysize(di):
    di = scprep.filter.filter_library_size(di, percentile=20, keep_cells='above')
    di = scprep.filter.filter_library_size(di, percentile=80, keep_cells='below')
    disease_filter = filter_rare_genes(di)
    return disease_filter

def filter_rare_genes(df):
    move_col = [ str(col) for col in df.columns if len(df[df[col] == 0][col]) >= len(df.index) * 0.5]
    df = df.drop(move_col,axis=1)
    return df

def filter_by_background(df,bg):
    columns  = df.columns
    columns_back = bg.columns
    p_col = []
    df_back = pd.DataFrame()
    for column in columns:
        if column in columns_back:
            wil,p_wil = stats.mannwhitneyu(df[column],bg[column],alternative='two-sided',use_continuity=False)
            p_col.append({"column":column,"pvalue":p_wil})
        else:
            p_col.append({"column":column,"pvalue":0.005})
    col_not = [i["column"] for i in p_col if i["pvalue"] >=0.5]
    df_back = df.drop(col_not,axis=1)
    return df_back

#新样本处理
def embedding(tab):
    tab = normal_data_processing(tab)
    tab = pre_embedding_data(tab)
    # 认为library size 小于1的部分没有意义
    #### 先根据自己的数据特征去掉一部分
    ##### 看motif是不是在大多数样本中都有表达,这里假定的是超过样本数目百分之50的样本没有表达
    ### filter rare motifs
    ##### 使用wilcoxon检验根据background删除一些特征
#### 找到这些不显著的特征AA，并且在相对应的疾病中删除
#### 假设服从同一分布，如果p值小于0.05，说明原假设不成立，则为显著性
    if tab.shape[0] == 1:
        pass
    else:
        tab = filter_by_librarysize(tab)
        tab = filter_rare_genes(tab)
        tab = filter_by_background(tab,background)
    tab.index.name = None
    tab.columns.name = None
    return tab

def kmer(tab):
    #补全所有的3kmer motif
    AAlist = "AFCUDNEQGHLIKOMPRSTVWY"
    Fmers = list(set([ f'{x}{y}{z}' for x in AAlist for y in AAlist for z in AAlist ]))
    val = list((np.zeros(tab.shape[0])))
    add_kmer = pd.DataFrame()
    for kmer in set(Fmers) - set(tab.columns):
        add_kmer[kmer] = val
    add_kmer.index = tab.index
    tab = pd.concat([tab,add_kmer],axis=1)
    tab.fillna(0, inplace=True)
    return tab

def intergration(df_li,df,name):
    df_li.append(df)
    name.append("New")
    exp, sample_labels = scprep.utils.combine_batches(
        df_li, 
        name,
        append_to_cell_names=True,
        common_columns_only=False)
    #消除深度影响
    exp = scprep.normalize.library_size_normalize(exp)
    #z-score normalize
    ind = exp.index
    col = exp.columns
    exp_scaled = pd.DataFrame(preprocessing.scale(exp),index=ind,columns=col)
    #补全kmer
    exp_scaled = kmer(exp_scaled)
    exp_scaled = exp_scaled.sort_index(axis=1)
    #只要输入样本的exp
    exp_scaled = exp_scaled.tail(1)
    return (exp_scaled,sample_labels)


@embedding_api.route("/emdeddingAnalysis",methods=['POST','GET'])
def ed():
    try:
        print ("embedding start")
        em_file_path = tempfile.NamedTemporaryFile(prefix="embedding", dir="/tmp").name
        embedding_file_path = (em_file_path + ".png")
        disease = [covid,nsclc,mel,breast,cmv,chy,cro,yel]
        # background = pd.read_hdf("/workspace2/yuet/tcrosetta/background_data/background_motif.hdf")
        name_list =["COVID-19", "NSCLC","Melnoma","Breast Cancer","CMV","classical Hodgkin lymphoma","Crohn’s disease","Yellow fever vaccine"]
        req = request.json
        fp = req["fp"]
        tab = pd.read_csv(fp,index_col=0)
        print ("start data progress")
        tab = embedding(tab)
        print ("embedding new data finish!")   
        tab,sample_labels = intergration(disease,tab,name_list)
        print ("embedding intergration finish")
        print ("embedding start")
        tab_phate = loaded_model.transform(tab)
        Y_umap = np.concatenate([load_phate,tab_phate])
        ss = scprep.plot.scatter2d(Y_umap,c = sample_labels,figsize=(14,8), cmap="Spectral",ticks=False,label_prefix="PHATE")
        plt.axhline(list(Y_umap[-1])[1],color="black",alpha=0.4,lw=2)
        plt.axvline(list(Y_umap[-1])[0],color="black",alpha=0.4,lw=2)
        plt.savefig(embedding_file_path)
        #就不计算中心点了，没必要，把图放上去给他们看看就算了
        print ("embedding save fig finish")
        with open(embedding_file_path, "rb") as handle:
            image_info = base64.b64encode(handle.read()).decode()
        
        os.system("rm -f %s" % embedding_file_path)
        result = {
            "embedding_image":image_info
        }
    except:
        result = {
            "status":400
        }
    #mycol.insert_one({"analysis":"embedding","data":result})
    print ("embedding finish")
    return json.dumps(result)
 