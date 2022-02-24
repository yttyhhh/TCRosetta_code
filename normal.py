
import os
import sys
import uuid
import json
import numpy as np
import pandas as pd
from bson.json_util import dumps
import base64
import matplotlib.pyplot as plt
import logomaker
from flask_cors import CORS
from flask import Flask,request,Blueprint
from Example import mycol

normal_api = Blueprint('normal', __name__)


# class MyEncoder(json.JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, bytes):
#             return str(obj,encoding='utf-8')
#         elif isinstance(obj,numpy.int64):
#             return int(obj)
# #         elif isinstance(obj, np.integer):
# #             return int(obj)
#         else:
#             return json.JSONEncoder.default(self, obj)

def normal_data_processing(df):
    #删除缺失数据
    df = df.dropna(subset=["AASeq","cloneFraction","Vregion","Dregion","Jregion"])
    #cloneFraction小于某个阈值的
    df = df[ df['cloneFraction'] > 0.00001]
    df["Condition"] = "New"
    if df.shape[0] == 0:
        df = pd.DataFrame()
    else:
    #CDR3格式筛选
        df = df[df["AASeq"].apply(lambda x: x[0] == "C" and x[-1] == "F" and "*" not in x and "_" not in x)]
        df = df[df["AASeq"].apply(lambda x: len(x) >= 8 and len(x) <=20)]
        #新加两列
        df['Count'] = df.groupby("AASeq")['cloneFraction'].transform('count')
        df['Ratio'] = df.groupby('AASeq')['cloneFraction'].transform('sum')
        #去重
        df = df.drop_duplicates(subset=['AASeq'],keep='first')
    return df


def logo(tab):
    figure,ax = plt.subplots(2,1, figsize=(10, 6))
    figure.subplots_adjust(wspace =0, hspace = 0.4)
    for iax, (_, pos, title) in enumerate([ ('front', [0, 1, 2, 3, 4], 'Position in CDR3 from start'), ('back', [-5, -4, -3, -2, -1], 'Position in CDR3 to end') ]):
        maker = []
        for idx, item in tab.iterrows():
            freq = item['cloneFraction']
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



@normal_api.route("/logo",methods=['POST'])
def common_analysis():
    # req = request.json
    # for item in req['file_list']:
    #     print(item)
    #     tab = pd.read_csv(item['response']["filePath"], index_col=0)
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    if tab.shape[0] == 0:
        tab = []
        return tab
    else:
        result = {
            "logo_image":logo(tab)
        }
        #mycol.insert_one({"analysis":"logo","data":result})
        return dumps(result)