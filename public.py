import sys
import math
import json
import pandas as pd
import tempfile
import os
import olga.load_model as load_model
import olga.generation_probability as pgen
import olga.sequence_generation as seq_gen
import scipy.stats as stats
from flask_cors import CORS
from flask import Flask,request, Blueprint
from pyecharts.charts import Boxplot
from jinja2 import Environment, FileSystemLoader
from pyecharts.globals import CurrentConfig
from pyecharts import options as opts
from normal import normal_data_processing
from Example import mycol


CurrentConfig.GLOBAL_ENV = Environment(loader=FileSystemLoader("./templates"))

public_api = Blueprint('public', __name__)

#这里的refernece用libo他们那个跑一下，看看还剩下多少序列，然后用这些序列来打分，作为这里的background分布
refernece_50 = pd.read_hdf("/workspace2/yuet/tcrosetta/background_data/background_prob_5.hdf")
params_file_name = '/workspace2/yuet/tools/OLGA/olga/default_models/human_T_beta/model_params.txt'
marginals_file_name = '/workspace2/yuet/tools/OLGA/olga/default_models/human_T_beta/model_marginals.txt'
V_anchor_pos_file ='/workspace2/yuet/tools/OLGA/olga/default_models/human_T_beta/V_gene_CDR3_anchors.csv'
J_anchor_pos_file = '/workspace2/yuet/tools/OLGA/olga/default_models/human_T_beta/J_gene_CDR3_anchors.csv'
#Load data
genomic_data = load_model.GenomicDataVDJ()
genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
#Load model
generative_model = load_model.GenerativeModelVDJ()
generative_model.load_and_process_igor_model(marginals_file_name)
#Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)


def public_analysis(df):
    # #Define the files for loading in generative model/data
    re = []
    prob = []
    if "Vregion" in list(df.columns):
        for idx,ins in df.iterrows():
            score = pgen_model.compute_aa_CDR3_pgen(ins["AASeq"],ins["Vregion"],ins["Jregion"])
            re.append({"AASeq":ins["AASeq"],"Vregion":ins["Vregion"],"Jregion":ins["Jregion"],"cloneFraction":ins["cloneFraction"],"Score":score})
            prob.append(score)
    else:
        for idx,ins in df.iterrows():
            score = pgen_model.compute_regex_CDR3_template_pgen(ins["AASeq"])
            re.append({"AASeq":ins["AASeq"],"Score":score})
            prob.append(score)
    re = pd.DataFrame(re)
    U1,p = stats.mannwhitneyu(refernece_50.sample(n=1000,axis=0)["Prob"],prob,alternative='two-sided',use_continuity=False)
    re["Pvalue"] = p
    score3 = [-math.log10(i) for i in prob if i != 0]
    back = []
    for i in refernece_50.Prob.values:
        if i == 0.0:
            pass
        else:
            back.append(-math.log10(i))
    box = Boxplot()
    data = box.prepare_data([score3])
    back_data = box.prepare_data([back])
    if "Vregion" in list(re.columns):
        tabData = [ {"AASeq":ins["AASeq"],"Vregion":ins["Vregion"],"Jregion":ins["Jregion"],"cloneFraction":'{:g}'.format(ins["cloneFraction"]),"Score":'{:g}'.format(ins["Score"])} for idx,ins in re.iterrows()]
        scdata = [ [-math.log10(ins["Score"]),-math.log10(ins["cloneFraction"])] for idx,ins in re.iterrows() if ins["cloneFraction"] and ins["Score"] != 0]
    else:
        tabData = [ {"AASeq":ins["AASeq"],"Score":'{:g}'.format(ins["Score"])} for idx,ins in re.iterrows() ]
        scdata = []
    result = {
        "chart":data,
        "chart_back":back_data,
        "Pvalue":p,
        "tabData":tabData,
        "scdata":scdata}

    return result

@public_api.route("/publicAnalysis",methods=['POST','GET'])
def pa():
    print ("public analysis start")
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    print ("public file read finish")
    if "Vregion" in list(tab.columns):
        tab = normal_data_processing(tab)
    else:
        tab = tab[tab["AASeq"].apply(lambda x: x[0] == "C" and x[-1] == "F" and "*" not in x and "_" not in x)]
        tab = tab[tab["AASeq"].apply(lambda x: len(x) >= 8 and len(x) <=20)]
        tab = tab.drop_duplicates(subset=["AASeq"],keep='first')
    if tab.shape[0] == 0:
        return []
    else:
    #这里再分析的时候直接卡死阈值，上线数据就是1000条，超过1000条就很慢了，或者只给network中结果node的public分析
        if tab.shape[0] <= 1000:
            new_sample = tab
        else:
            new_sample = tab.sample(n=1000,random_state=None,replace=False,axis=0)
        try:
            pub = public_analysis(new_sample)
        except:
            pub = {
                "status":400
            }
        print ("public analysis finish")
        #mycol.insert_one({"analysis":"public","data":pub})
        return json.dumps(pub)