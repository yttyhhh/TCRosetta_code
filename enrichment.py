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
from fisher import pvalue
from Example import mycol

enrichment_api = Blueprint('enrichment', __name__)

# global enrichment_file_path1
sequence_count = pd.read_csv("/workspace2/yuet/tcrosetta/sequence_count.csv",sep=",",index_col=0)
small_size = pd.read_csv("/workspace2/yuet/tcrosetta/enrichment/disease_with_low_sequence.csv")

#这个是用来去画pie的。相对数量
def relative_number(es_result,back=sequence_count):
    tab = pd.DataFrame(es_result)
    relative_number = []
    for cond,number in Counter(tab["Condition"]).most_common()[0:10]:
        background_number = int(back[back["Condition"] == cond]["sequencenumber"].values)
        zb = '{:g}'.format((number / background_number))
        relative_number.append({"value":zb,"name":cond})
    relative_number = sorted(relative_number,key=lambda x:x["value"],reverse=True)
    return relative_number

def enrichment_cluster(tab,filepath):
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
    os.system("rm %s" % filepath)
    os.system("rm %s" % open_dir)
    return clus

#模糊匹配，小于1000条时采用
def esearch1000(sequence):
    es = Elasticsearch(
        '127.0.0.1:9200',
        sniff_on_start=True,
        sniff_on_connection_fail=True,
        sniffer_timeout=60)
    if es.ping():
        CDR_list = []
        AASeq_contain = {
            "from":0,
            "size":1000,
            'query':{
                "fuzzy":{
                    "AASeq":{
                        "value":sequence.lower(),
                        "prefix_length":2,
                        "fuzziness":1
                    }
                }
            }
        }
        result = es.search(
            index="elasticsearch",
            doc_type="es",
            body=AASeq_contain,
            request_timeout = 60)
        if result['hits']['hits'] == []:
            print ("no result")
            return CDR_list
        else:
            for seq in result['hits']['hits']:
                CDR_list.append(seq['_source'])
            return  CDR_list
    else:
        print ("error in annnotation")
        CDR_list = []
        return CDR_list

#精确匹配大于1000 小于3000时使用
def esearch3000(sequence):
    es = Elasticsearch(
        '127.0.0.1:9200',
        sniff_on_start=True,
        sniff_on_connection_fail=True,
        sniffer_timeout=60)
    if es.ping():
        CDR_list = []
        AASeq_contain = {
            "from":0,
            "size":1000,
            'query':{
                "match":{
                    "AASeq":sequence.lower()
                }
            }
        }
        result = es.search(
            index="elasticsearch",
            doc_type="es",
            body=AASeq_contain,
            request_timeout = 60)
        if result['hits']['hits'] == []:
            print ("no result")
            return CDR_list
        else:
            for seq in result['hits']['hits']:
                CDR_list.append(seq['_source'])
            return  CDR_list
    else:
        print ("error in annnotation")
        CDR_list = []
        return CDR_list

def enrichment(es_result):
    if len(es_result) == []:
        final = []
        # final.to_csv(filepath)
        return final
    else:
        final = pd.DataFrame(es_result)
        #table data
        final = final[["AASeq","Vregion","Dregion","Jregion","Condition"]]
        # final.to_csv(filepath)
        condition = [ i["Condition"] for i in es_result ]
        count_in_sample = Counter(condition)
        count_in_sample = count_in_sample.most_common()
        if len(count_in_sample) > 20:
            dat = [{"name":str(idx),"children":[],"value":int(ins)} for idx,ins in count_in_sample[0:20]]
        else:
            dat = [{"name":str(idx),"children":[],"value":int(ins)} for idx,ins in count_in_sample]
        for cond in dat:
            cr = final[final["Condition"] == cond["name"]]
            count_cr = Counter(cr["Vregion"]).most_common()
            if len(count_cr) > 20:
                chil_v = [{"name":key,"value":value} for key,value in count_cr[0:20]]
            else:
                chil_v = [{"name":key,"value":value} for key,value in count_cr]
            v = {
                "name":"Vregion",
                "value":6,
                "children":chil_v}
            count_jr = Counter(cr["Jregion"]).most_common()
            if len(count_jr) > 20:
                chil_j = [{"name":key,"value":value} for key,value in count_jr[0:20]]
            else:
                chil_j = [{"name":key,"value":value} for key,value in count_jr]
            j = {
                "name":"Jregion",
                "value":4,
                "children":chil_j}
            count_cdr3 = Counter(cr["AASeq"]).most_common()
            if len(count_cdr3) > 50:
                chil_cdr3 = [{"name":key,"value":value} for key,value in count_cdr3[0:50]]
            else:
                chil_cdr3 = [{"name":key,"value":value} for key,value in count_cdr3]
            cdr3 = {
                "name":"CDR3",
                "value":10,
                "children":chil_cdr3}
            cond['children'] = [v,j,cdr3]
        tabData = final
        re = [ {"AASeq":ins["AASeq"],"Vregion":ins["Vregion"],"Dregion":ins["Dregion"],"Jregion":ins["Jregion"],"Condition":ins["Condition"]} for idx,ins in tabData.iterrows()]
        # result = {
        #     "chart":dat,
        #     "tabData":re
        # }
        return dat , re

def specificity(df,ss=small_size):
    df_re = []
    if df.shape[0] == 2:
        if "COVID-19" not in list(df["Condition"].values):
            df = df.sort_values(by="Count",ascending=False)
            z = int(df.iloc[[0]]["Count"])
            f = int(df.iloc[[1]]["Count"])
            #这里要加一个样本数目少的疾病的list，判断其是否在该list中
            if z >= f * 2:
                df_re.append(df.iloc[[0]])
                if df.iloc[[1]]["Condition"].values[0] in ss:
                    df_re.append(df.iloc[[1]])
                else:
                    pass
            else:
                df_re.append(df)
        else:
            df = df.sort_values(by="Count",ascending=False)
            z = int(df.iloc[[0]]["Count"])
            f = int(df.iloc[[1]]["Count"])
            if str(df.iloc[[0]]["Condition"]) == "COVID-19":
                if z >= f * 6 and int(df.iloc[[0]]["Count"]) > 30:
                    df_re.append(df.iloc[[0]])
                else:
                    df_re.append(df.iloc[[1]])
            else:
                if z >= f * 2:
                    df_re.append(df.iloc[[0]])
                else:
                    df_re.append(df)
    elif df.shape[0] > 2:
        #直接保留前两个，如果有COVID存在的话还是要先判断下其是否大于6个
        if "COVID-19" in list(df["Condition"].values):
            #大于6的话保留COVID-19,小于的话就可以掠过
            df = df.sort_values(by="Count",ascending=False)
            #如果covid在其中，并且coivd的数量大于其他加起来，则定位该序列就是covid的序列
            sum_no_covid = sum(df[df["Condition"] != "COVID-19"]["Count"])
            count_of_covid = int(df[df["Condition"] == "COVID-19"]["Count"])
            if count_of_covid >= sum_no_covid and count_of_covid > 30:
                co = df[df["Condition"] == "COVID-19"]
                df_re.append(co)
            else:
                co_el = df[df["Condition"] != "COVID-19"]
                co_el = co_el.sort_values(by="Count",ascending=False)
                co_first = int(co_el.iloc[[0]]["Count"])
                co_two = int(co_el.iloc[[1]]["Count"])
                cooo = co_el.iloc[[0]]
                df_re.append(cooo)
                if co_first > co_two * 2 and str(co_el.iloc[[1]]["Condition"]) not in ss:
                    pass
                else:
                    cooo_two = co_el.iloc[[1]]
                    df_re.append(cooo_two)
        else:
            #如果COVID-19不在其中的话，先排序
            df = df.sort_values(by="Count",ascending=False)
            if df.shape[0] == 3:
                if int(df.iloc[[0]]["Count"]) == int(df.iloc[[1]]["Count"])  == int(df.iloc[[2]]["Count"]) :
                    df_re.append(df)
                else:
                    coo = df.iloc[[0]]
                    f_first = int(df.iloc[[0]]["Count"])
                    #对于排在第二的，如果差别小于2倍，就留下，如果大于两倍就不要了
                    f_two = int(df.iloc[[1]]["Count"])
                    f_third = int(df.iloc[[2]]["Count"])
                    df_re.append(coo)
                    if f_first > f_two * 2 and str(df.iloc[[1]]["Condition"]) not in ss:
                        pass
                    else:
                        coo_two = df.iloc[[1]]
                        df_re.append(coo_two)
                        if f_two > f_third * 2 and str(df.iloc[[2]]["Condition"]) not in ss:
                            pass
                        else:
                            coo_third = df.iloc[[2]]
                            df_re.append(coo_third)
            #取排在第一的
            else:
                coo = df.iloc[[0]]
                f_first = int(df.iloc[[0]]["Count"])
                #对于排在第二的，如果差别小于2倍，就留下，如果大于两倍就不要了
                f_two = int(df.iloc[[1]]["Count"])
                df_re.append(coo)
                if f_first > f_two * 3 and str(df.iloc[[1]]["Condition"]) not in ss:
                    pass
                else:
                    coo_two = df.iloc[[1]]
                    df_re.append(coo_two)
    return df_re

def data_for_specificity(cdr3_list,ss=small_size):
    ss = list(ss["Condition"].values)
    final = pd.DataFrame(cdr3_list)
    #table data
    final = final[["AASeq","Condition"]]
    final["Status"] = "Unknown"
    final["Count"] = final.groupby(["AASeq","Condition"])["Status"].transform("count")
    etra_cond = []
    for idx,ins in final.iterrows():
        if ins["Count"] == 1 and ins["Condition"] in ss:
            etra_cond.append(ins)
    etra_cond1 = pd.DataFrame(etra_cond)
    final = final[final["Count"] > 1]
    final = pd.concat([final,etra_cond1],axis=0)
    final = final.sort_values(by=["AASeq","Condition"])
    final = final.drop_duplicates(subset=['AASeq',"Condition"],keep='first')
    final = final[final["Condition"] != "Healthy"]
    seq = set(list(final["AASeq"]))
    ww = []
    pp = []
    for ii in seq:
        dd = final[final["AASeq"] == ii]
        if dd.shape[0] > 3:
            pass
        elif dd.shape[0] == 1:
            pp.append(dd)
        else:
            ww.append(dd)
    ff_first = pd.concat(ww,axis=0)
    ff_seq = set(list(ff_first["AASeq"]))
    result = []
    for i in ff_seq:
        subdf = ff_first[ff_first["AASeq"] == i]
        co_result = specificity(subdf,small_size)
        result.extend(co_result)
    multi_seq_result = pd.concat(result,axis=0)
    mm = pd.concat(pp,axis=0)
    mm_no_covid = mm[(mm["Condition"] != "COVID-19")&(mm["Condition"] != "Non-small cell lung cancer")&(mm["Condition"] != "T1D")]
    covid_filter = mm[(mm["Condition"] == "COVID-19")&(mm["Count"] > 30)]
    nsclc_filter = mm[(mm["Condition"] == "Non-small cell lung cancer")&(mm["Count"] > 3)]
    t1d_filter = mm[(mm["Condition"] == "T1D")&(mm["Count"] > 3)]
    ff_two = pd.concat([mm_no_covid,covid_filter,nsclc_filter,t1d_filter],axis=0)
    first_two_result = pd.concat([multi_seq_result,ff_two])
    result = []
    for idx,ins in Counter(first_two_result["Condition"]).most_common()[0:10]:
        result.append([str(idx),int(ins)])
    specific_xdata = []
    specific_ydata = []
    color_plat = ["#1C0C5B","#142850","#27496D","#32407B","#3D2C8D","#664E88","#0C7B93","#63B4B8","#94D3AC","#CCEDD2"]
    ind = 0
    for idx,ins in result:
        specific_xdata.append(idx)
        # if ins == result[0][1]:
        #     specific_ydata.append({"value":ins,"itemStyle":{"color":'#a90000'}})
        # else:
        specific_ydata.append({"value":ins,"itemStyle":{"color":color_plat[ind]}})
        ind = ind + 1
    # return {
    #     "funnel":result,
    #     "fun_legend":legend_result
    # }
    result = {
        "specificity_xdata":specific_xdata,
        "specificity_ydata":specific_ydata
    }
    return result

def p_test(finalr,sc = sequence_count):
    sc = sc.sort_values(by=["sequencenumber"],ascending=False)
    final_df = pd.DataFrame(finalr)
    final_df = final_df[["AASeq","RunId","Condition"]]
    final_df = final_df.drop_duplicates(subset=['AASeq',"RunId","Condition"],keep='first')
    final_df["Status"] = "Unknown"
    final_df["Count"] = final_df.groupby(["AASeq","Condition"])["Status"].transform("count")
    final_df = final_df.drop_duplicates(subset=["AASeq","Condition"],keep='first')
    condition = [ ins["Condition"] for idx,ins in final_df.iterrows() ]
    condition = set(condition)
    N = 244955895
    n = sum(final_df["Count"])
    presult = []
    condition_bg = [str(i) for i in sc["Condition"]]
    for i in condition:
        if i in condition_bg:
            M = int(sc[sc["Condition"] == i]["sequencenumber"].values[0])
            k = sum(final_df[final_df["Condition"] == i]["Count"])
            p = pvalue(M-k, N-M-n+k, k, n-k)
            pval = {"Condition":i,"Pvalue":p.left_tail}
            presult.append(pval)
        else:
            pass
    return presult
    
    

#network enrichment
@enrichment_api.route("/enrichment",methods=['POST','GET'])
def en():
    try:
        print ("enrichment start")
        enrichment_file_path = tempfile.NamedTemporaryFile(prefix="enrichment", dir="/tmp").name
        # if request.method == 'GET':
        #这部分做enrichment的需要在前端先判断下network完成了，才能进行下一步
        #这里的路径需要考虑到network
        req = request.json
        fp =req["fp"]
        tab = pd.read_csv(fp)
        if "Vregion" in list(tab.columns):
            tab = normal_data_processing(tab)
            tab = enrichment_cluster(tab,enrichment_file_path)
            print ("enrichment cluster finish")
        else:
            tab = tab[tab["AASeq"].apply(lambda x: x[0] == "C" and x[-1] == "F" and "*" not in x and "_" not in x)]
            tab = tab[tab["AASeq"].apply(lambda x: len(x) >= 8 and len(x) <=20)]
            tab = tab.drop_duplicates(subset=["AASeq"],keep='first')
        if tab.shape[0] > 1000 and tab.shape[0] <= 3000:
            # tab = tab.sample(n=1000,random_state=None,replace=False,axis=0)
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
        data_for_pie = relative_number(final_result)
        print ("elasticsearch search finish")
        enrichment_chart,enrichment_data = enrichment(final_result)
        select_data = list(set([ i["Condition"] for i in enrichment_data]))
        select_data = [ {"value": i,"label":i} for i in select_data]
        print ("enrichment finish")
        print ("specificity start")
        specificity_p = p_test(final_result,sequence_count)
        # specificity_data= data_for_specificity(cdr3_list=final_result,ss=small_size)
        # pr = [name for name in specificity_p if name["Condition"] in specificity_data["specificity_xdata"]]
        pr = [name for name in specificity_p]
        pr.sort(key=lambda k :k["Pvalue"],reverse=False)
        pr = [{"Condition":i["Condition"],"Pvalue":'{:g}'.format(i["Pvalue"])} for i in pr if float(i["Pvalue"]) < 0.05]
        result = {
            "enrichment_data":enrichment_data,
            "select_data":select_data,
            "enrichment_chart":enrichment_chart,
            # "specificity_data":specificity_data,
            "specificity_p":pr,
            "relative_number":data_for_pie
        }
    except:
        result = {
            "status": 400
        }
    #mycol.insert({"analysis":"enrichment","data":result})
    print ("specificity finish")
    return json.dumps(result)


    
