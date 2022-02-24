from flask_cors.core import set_cors_headers
import pandas as pd
import json
import os
import tempfile
from collections import Counter
from elasticsearch import Elasticsearch
from operator import itemgetter
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import logo, normal_data_processing
from enrichment import enrichment_cluster,esearch1000,esearch3000
from Example import mycol

onlysearch_api = Blueprint('onlysearch', __name__)

@onlysearch_api.route("/onlysearch",methods=['POST'])
def only():
    print ("only search and annotation start")
    only_search_file_path = tempfile.NamedTemporaryFile(prefix="onlySearch", dir="/tmp").name
    req = request.json
    fp =req["fp"]
    tab = pd.read_csv(fp)
    #判断有无Vregion来判断是输入序列还是上传文件
    if "Vregion" in list(tab.columns):
        tab = normal_data_processing(tab)
        tab = enrichment_cluster(tab,only_search_file_path)
        print ("only search cluster by GIANA finish")
    else:
        tab = tab[tab["AASeq"].apply(lambda x: x[0] == "C" and x[-1] == "F" and "*" not in x and "_" not in x)]
        tab = tab[tab["AASeq"].apply(lambda x: len(x) >= 8 and len(x) <=20)]
        tab = tab.drop_duplicates(subset=["AASeq"],keep='first')
    if tab.shape[0] > 1000 and tab.shape[0] <= 3000:
        sample_list = list(tab["AASeq"])
        sample_list = [str(seq.lower()) for seq in sample_list]
        final_result = []
        for seq in sample_list:
            final_result.extend(esearch3000(seq))
    elif tab.shape[0] <= 1000:
        sample_list = list(tab["AASeq"])
        sample_list = [str(seq.lower()) for seq in sample_list]
        final_result = []
        for seq in sample_list:
            final_result.extend(esearch1000(seq))
    else:
        tab = tab.sample(n=3000,random_state=None,replace=False,axis=0)
        sample_list = list(tab["AASeq"])
        sample_list = [str(seq.lower()) for seq in sample_list]
        final_result = []
        for seq in sample_list:
            final_result.extend(esearch3000(seq))
    final = pd.DataFrame(final_result)
    final = final[["AASeq","Vregion","Dregion","Jregion","Condition"]]
    re = [ {"AASeq":ins["AASeq"],"Vregion":ins["Vregion"],"Dregion":ins["Dregion"],"Jregion":ins["Jregion"],"Condition":ins["Condition"]} for idx,ins in final.iterrows()]
    
    select_data = list(set([ i["Condition"] for i in re]))
    select_data = [ {"value": i,"label":i} for i in select_data]
    print ("only search by elasticsearch finish number: %s" % len(final_result))
    result = {
        "only_search_data":re,
        "only_select_data":select_data
    }
    print ("only search finish")
    #mycol.insert_one({"analysis":"onlysearch","data":result})
    return json.dumps(result)

