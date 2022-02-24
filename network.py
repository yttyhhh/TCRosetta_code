from inspect import ClosureVars
import math
import os
import sys
import json
import tempfile
import pandas as pd
import numpy as np
from collections import Counter
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
# from pyecharts import options as opts
# from pyecharts.charts import Graph
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import normal_data_processing
import Levenshtein
import base64
import uuid
import matplotlib.pyplot as plt
import logomaker
from Example import mycol

network_api = Blueprint('network', __name__)
#1.normal data processing
#2.创建可以作为输入GIANA的的文件 entwork.tsv
#3.运行GIANA得到结果文件entwork--RotationEncodingBL62.txt
#4.提取结果文件中的序列，并进行比对 pre_alignment.fasta
#5.读取比对后的结果文件，计算编辑距离 post_alignment.fasta
#6.保存距离矩阵，并且作为R脚本的输入
#7.得到R脚本的输出node.csv,link.csv，
#8.画图
def network_logo(tab):
    figure,ax = plt.subplots(2,1, figsize=(10, 6))
    figure.subplots_adjust(wspace =0, hspace = 0.4)
    for iax, (_, pos, title) in enumerate([ ('front', [0, 1, 2, 3, 4, 5, 6, 7], 'Position in CDR3 from start'), ('back', [-7,-6,-5, -4, -3, -2, -1], 'Position in CDR3 to end') ]):
        maker = []
        for idx, item in tab.iterrows():
            freq = item['Ratio']
            seq = item['AASeq']
            for idx in pos:
                maker.append( [idx+1, seq[idx], freq] )
        maker = pd.DataFrame(maker, columns=['Pos', 'Chr', 'Frequency'])
        make = maker.pivot_table(index='Pos', columns='Chr', values='Frequency')
        make.fillna(0, inplace=True)
        for i in make.index:
            make.loc[i, :] = make.loc[i, :] / sum(make.loc[i, :])
        #print(make)
        ax[iax].set_xlabel(title)
        ww_logo = logomaker.Logo(make,
                         font_name='Stencil Std',
                         color_scheme='NajafabadiEtAl2017',
                         ax = ax[iax],
                         vpad=.1,
                         width=1)
    temp_name = uuid.uuid1()
    figure.savefig(f"{temp_name}.png", bbox_inches='tight')
    figure.clf()
    plt.close()
    with open(f"{temp_name}.png", "rb") as handle:
        res = base64.b64encode(handle.read()).decode()
    os.system(f"rm -f {temp_name}.png")
    return res

def cluster(tab,filepath):
    # save data to do GIANA
    d = tab.pop("Vregion")
    tab.insert(1,"Vregion",d)
    tab = tab[["AASeq","Vregion","Jregion","Ratio"]].rename(columns={"AASeq":"aminoAcid","Vregion":"vMaxResolved"})
    tab.to_csv(filepath,sep="\t",index=0)
    order = ("python /workspace2/yuet/tcrosetta/GIANA/GIANA4.py -f %s -o /tmp/" % filepath)
    os.system(order)
    #read GIANA result
    aaseq = []
    cluster_id = []
    vregion = []
    jregion = []
    ratio = []
    open_dir = ("/tmp/%s--RotationEncodingBL62.txt" % filepath.split("/")[2])
    with open(open_dir,"r") as handle:
        for line in handle:
            if "#" in line:
                pass
            else:
                aaseq.append(line.split("\t")[0])
                cluster_id.append(line.split("\t")[1])
                vregion.append(line.split("\t")[2])
                jregion.append(line.split("\t")[3])
                ratio.append(float(line.split("\t")[4]))
    clus = pd.DataFrame()
    clus["AASeq"] = aaseq
    clus["cluster_id"] = cluster_id
    clus["Vregion"] = vregion
    clus["Jregion"] = jregion
    clus["Ratio"] = ratio
    
    #network logo(cluster result logo)
    # network_logo_base = network_logo(clus)

    # save GIANA result to display in table
    # clus.to_csv(cluster_file_path,index=0)
    # make fastq format data to do alignment
    alignment_file_path = (filepath + "_" + "alignment.fasta")
    with open(alignment_file_path,"w") as handle:
        for idx,ins in clus.iterrows():
            title = (">sequence/" + str(idx))
            handle.write(title)
            handle.write("\n")
            handle.write(str(ins["AASeq"]))
            handle.write("\n")
    #do alignment with muscle
    alignment_out_file_path = (filepath + "_" + "alignment.out.fasta")
    muscle_order = ("/data/yuet/R/tool/muscle3.8.31_i86linux64 -in %s -out %s" % (alignment_file_path,alignment_out_file_path))
    os.system(muscle_order)
    # read alignment result
    ind = []
    aa = []
    with open(alignment_out_file_path,"r") as handle:
        for line in handle:
            if ">" in line:
                idx = int(line.split("/")[1].split("\n")[0])
                ind.append(idx)
            else:
                seq = str(line.split("\n")[0])
                aa.append(seq)
    result = pd.DataFrame()
    result["Index"] = ind
    result["AASeq"] = aa
    result = result.sort_values(by=["Index"],ascending=True)
    mat_dis = []
    for i in result.AASeq.values:
        distance = []
        for j in result.AASeq.values:
            ed = Levenshtein.distance(i,j)
            if ed <= 3:
                distance.append(1)
            elif ed == 0:
                distance.append(0)
            else:
                distance.append(0)
        mat_dis.append(distance)
    mat_dis = np.array(mat_dis)
    row,col = np.diag_indices_from(mat_dis)
    mat_dis[row,col] = 0
    df_mat_dis = pd.DataFrame(mat_dis)
    df_mat_dis.columns = clus.AASeq.values
    df_mat_dis.index = clus.AASeq.values
    cluster_matrix_file_path = (filepath + "_"+"cluster" + "_" + "matrix.csv")
    df_mat_dis.to_csv(cluster_matrix_file_path,index=0)
    node_file_path = (filepath + "_" + "node.csv")
    link_file_path = (filepath + "_" + "link.csv")
    os.system("/home/yuet/anaconda3/envs/r410/bin/Rscript /data/yuet/R/script/get_node_edge.R %s %s %s" % (cluster_matrix_file_path,node_file_path,link_file_path))
    print ("igraph finish")
    node = pd.read_csv(node_file_path,index_col=0)
    link = pd.read_csv(link_file_path,index_col=0)
    nodes = [{"name":str(ins["Node"]),"symbolSize":300* float(ins["Weight"]),"itemStyle": {"normal": {"color": ins["Color"]}},} for idx,ins in node.iterrows()]
    links = [{"source":str(ins["Source"]),"target":str(ins["Target"])} for idx,ins in link.iterrows()]
    tabData = clus
    node = pd.read_csv(node_file_path,index_col=0)
    node =node.rename(columns={"Node":"AASeq"})
    ss = pd.merge(tabData,node,on="AASeq")
    re = [ {"AASeq":ins["AASeq"],"Vregion":ins["Vregion"],"Jregion":ins["Jregion"],"Ratio":'{:g}'.format(ins["Ratio"]),"Weight":'{:g}'.format(ins["Weight"]),"Cluster":ins["Cluster"]} for idx,ins in ss.iterrows()]
    re = sorted(re,key=lambda x : x["Cluster"],reverse=False)
    print ("start draw logo")
    #network cluster different logo
    cluster_id = set(ss["Cluster"])
    print (cluster_id)
    cluster_of_sample = [ ss[ss["Cluster"] == i] for i in cluster_id if ss[ss["Cluster"] == i].shape[0] > 5 ]
    # print (len(cluster_of_sample))
    all_network_base = [network_logo(i) for i in cluster_of_sample]
    # print (len(all_network_base))
    # print (all_network_base[0])
    network_logo_base = [{"value": ins,"label":"cluster"+ str(idx+1)} for idx,ins in enumerate(all_network_base)] 
    print (len(network_logo_base))
    # print ("enrichment finish")
    os.system("rm %s" % filepath)
    os.system("rm %s" % open_dir)
    os.system("rm %s" % node_file_path)
    os.system("rm %s" % link_file_path)
    os.system("rm %s" % cluster_matrix_file_path)

    result = {
        "nodes":nodes,
        "links":links,
        "tabData":re,
        "logoE":network_logo_base,
    }
    return result

@network_api.route("/network",methods=['POST','GET'])
def nt():
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    if tab.shape[0] == 0:
        tab = []
        return tab
    else:
        network_file_path = tempfile.NamedTemporaryFile(prefix="network", dir="/tmp").name
        try:
            result = cluster(tab,network_file_path)
        except:
            result = {
                "status" : 400
            }
        # mycol.insert_one({"analysis":"network","data":result})
    return json.dumps(result)