import os
import json
import pandas as pd
import re
from collections import Counter


def most_freq(item, freq):
    dd = Counter()
    for l, f in zip(item, freq):
        dd[l] = dd[l] + f
    return dd.most_common()[0][0]


def region_rename(x):
    """
    format V/J
    """
    if 'TRBV' in x:
        try:
            return "".join(re.findall(r"(TRBV[0-9]{1,2}-[0-9]{1,2})|(TRBV[0-9]{1,2})", x)[0])
        except:
            return x
    else:
        try:
            return "".join(re.findall(r"(TRBJ[0-9]{1,2}-[0-9]{1,2})|(TRBJ[0-9]{1,2})", x)[0])
        except:
            return x

def parser_data(file_path, data_type):
    """
    用来解析用户上传的TCR数据
    Input:
        file_path: [String], 输入文件上传之后的位置
        data_type: [String], 输入文件的数据格式，只能是 ["MiXCR", "IMSEQ", "RTCR", "CATT"] 中的一种

    Output:
        dataframe 格式，包含有两列 AASeq, cloneFraction
    """

    if data_type not in ["MiXCR", "IMSEQ", "RTCR", "CATT","TCRosetta","CDR3"]:
        print("Data type error when PARSER_DATA")

    res = []

    if data_type == "MiXCR":
        tab = pd.read_table(file_path, usecols=['aaSeqCDR3', 'cloneFraction', 'allJHitsWithScore', 'allVHitsWithScore'])
        tab = tab[ tab['aaSeqCDR3'].apply(lambda x: x[0]=='C' and x[-1]=='F' and "*" not in x and '_' not in x) ]
        tab['Vregion'] = tab['allVHitsWithScore'].apply(lambda x: x.split(',')[0].split("*")[0])
        tab['Jregion'] = tab['allJHitsWithScore'].apply(lambda x: x.split(',')[0].split("*")[0])

        for seq, group in tab.groupby("aaSeqCDR3"):
            res.append({
                          "AASeq": seq,
                          "cloneFraction": sum(group['cloneFraction']),
                          "Vregion": most_freq(group['Vregion'], group['cloneFraction']),
                          "Jregion": most_freq(group['Jregion'], group['cloneFraction'])
                        })

        res = pd.DataFrame(res)
        res = res.fillna('None')
        res = res[res["Vregion"] != "None"]
        res = res[res["Jregion"] != "None"]
        res = res[res["AASeq"] != "None"]
        res = res[res["cloneFraction"] != "None"]
        total = sum(res['cloneFraction'])
        res['cloneFraction'] = res['cloneFraction'] / total
        res["Dregion"] = "Unknown"
        # print (res.head())


    #TODO(chensy) Remind user to use correct command parameters
    if data_type == "IMSEQ":
        tab = pd.read_table(file_path, header=None)
        tab.columns = ['Org', 'cloneFraction']
        tab['Vregion'] = tab['Org'].apply(lambda x: 'TRBV' + x.split(':')[0].split("V")[1])
        tab['Jregion'] = tab['Org'].apply(lambda x: 'TRBJ' + x.split(':')[2].split("J")[1])
        tab['AASeq'] = tab['Org'].apply(lambda x: x.split(':')[1] )

        for seq, group in tab.groupby("AASeq"):
            res.append({
               "AASeq": seq,
               "cloneFraction": sum(group['cloneFraction']),
               "Vregion": most_freq(group['Vregion'], group['cloneFraction']),
               "Jregion": most_freq(group['Jregion'], group['cloneFraction'])
             })
        res = pd.DataFrame(res)
        res = res.fillna('None')
        res = res[res["Vregion"] != "None"]
        res = res[res["Jregion"] != "None"]
        res = res[res["AASeq"] != "None"]
        res = res[res["cloneFraction"] != "None"]
        total = sum(res['cloneFraction'])
        res['cloneFraction'] = res['cloneFraction'] / total
        res["Dregion"] = "Unknown"

    if data_type == "RTCR":
        tab = pd.read_table(file_path, usecols=['Amino acid sequence','Number of reads', 'V gene', 'J gene'])
        tab['Vregion'] = tab['V gene'].apply(lambda x: x.split('*')[0])
        tab['Jregion'] = tab['J gene'].apply(lambda x: x.split('*')[0])

        for seq, group in tab.groupby("AASeq"):
            res.append({
                "AASeq": seq,
                "cloneFraction": sum(group['cloneFraction']),
                "Vregion": most_freq(group['Vregion'], group['cloneFraction']),
                "Jregion": most_freq(group['Jregion'], group['cloneFraction'])
                })

        res = pd.DataFrame(res)
        res = res.fillna('None')
        res = res[res["Vregion"] != "None"]
        res = res[res["Jregion"] != "None"]
        res = res[res["AASeq"] != "None"]
        res = res[res["cloneFraction"] != "None"]
        total = sum(res['cloneFraction'])
        res['cloneFraction'] = res['cloneFraction'] / total
        res["Dregion"] = "Unknown"

    if data_type == "CATT":
        tab = pd.read_csv(file_path, usecols=['AAseq', 'Frequency', 'Vregion', 'Jregion'])
        tab['Vregion'] = tab['Vregion'].apply(lambda x: x.split('*')[0])
        tab['Jregion'] = tab['Jregion'].apply(lambda x: x.split('*')[0])
        for seq, group in tab.groupby("AAseq"):
            res.append({
                "AASeq": seq,
                "cloneFraction": sum(group['Frequency']),
                "Vregion": most_freq(group['Vregion'], group['Frequency']),
                "Jregion": most_freq(group['Jregion'], group['Frequency'])
                })

        res = pd.DataFrame(res)
        #fill NaN
        res = res.fillna('None')
        res = res[res["Vregion"] != "None"]
        res = res[res["Jregion"] != "None"]
        res = res[res["AASeq"] != "None"]
        res = res[res["cloneFraction"] != "None"]
        total = sum(res['cloneFraction'])
        res['cloneFraction'] = res['cloneFraction'] / total
        res["Dregion"] = "Unknown"

    #add by yuet
    if data_type == "TCRosetta":
        tab = pd.read_csv(file_path, usecols=['AAseq', 'Frequency', 'Vregion', 'Jregion'])
        tab['Vregion'] = tab['Vregion'].apply(lambda x: x.split('*')[0])
        tab['Jregion'] = tab['Jregion'].apply(lambda x: x.split('*')[0])
        for seq, group in tab.groupby("AAseq"):
            res.append({
                "AASeq": seq,
                "cloneFraction": sum(group['Frequency']),
                "Vregion": most_freq(group['Vregion'], group['Frequency']),
                "Jregion": most_freq(group['Jregion'], group['Frequency'])
                })

        res = pd.DataFrame(res)
        
        #fill NaN
        res = res.fillna('None')
        res = res[res["Vregion"] != "None"]
        res = res[res["Jregion"] != "None"]
        res = res[res["AASeq"] != "None"]
        res = res[res["cloneFraction"] != "None"]
        total = sum(res['cloneFraction'])
        res['cloneFraction'] = res['cloneFraction'] / total
        res["Dregion"] = "Unknown"

    #add by yuet
    if data_type == "CDR3":
        with open("/workspace/yuet/tcr_tool/example/EXAMPLE.CDR3","r") as handle:
            sequence = (str(handle.read()))
        acid = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
        status = 0
        res = pd.DataFrame()
        for seq in sequence:
            if seq == "\n":
                pass
            elif seq.upper() not in acid:
                status = 1
        if status == 0:
            rel_seq = []
            ss = sequence.split("\n")
            rel_seq = [ i for i in ss if len(i) < 20 and len(i) > 7 ]
        res["AASeq"] = rel_seq
    
    if "cloneFraction" in list(res.columns):
        res.sort_values('cloneFraction', ascending=True, inplace=True)
    return res
