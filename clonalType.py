import pandas as pd
import numpy as np
from jinja2 import Environment, FileSystemLoader
from pyecharts.globals import CurrentConfig
from pyecharts import options as opts
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import normal_data_processing
from pyecharts.charts import Pie
from operator import itemgetter
#from Example import mycol

CurrentConfig.GLOBAL_ENV = Environment(loader=FileSystemLoader("./templates"))

clonalType_api = Blueprint('clonalType', __name__)

def make_clonality(df):
    data_for_pie = []
    df = df.sort_values(by=["Ratio"],ascending=False)
    df_ten = df.head(10)
    total = sum(df["Ratio"])
    total_ten = sum(df_ten["Ratio"])
    other = total-total_ten
    # lev = np.percentile(np.array(df["Ratio"]),99.9)
    for idx,ins in df_ten.iterrows():
        # if ins["Ratio"] < lev:
        #     other = other + ins["Ratio"]
        # else:
        data_for_pie.append({"AASeq":str(ins["AASeq"]),"Ratio":ins["Ratio"]})
    data_for_pie.append({"AASeq":"Other","Ratio":other})
    data_for_pie = sorted(data_for_pie,key=itemgetter('Ratio'))
    data_for_pie = [ [i["AASeq"],i["Ratio"]] for i in data_for_pie]
    pie = (
        Pie()
        .add(
            series_name = "", 
            data_pair = data_for_pie,
            radius=["0%","80%"])
        .set_global_opts(
            title_opts=opts.TitleOpts(""),
            legend_opts=opts.LegendOpts(is_show=False),
            toolbox_opts = opts.ToolboxOpts(
                is_show=True,
                pos_top="3%",
                pos_right="right",
                feature={
                    "saveAsImage":{
                        "title":"save as image"
                    }
                }))
        .set_series_opts(
            label_opts = opts.LabelOpts(is_show=False),
            # tooltip_opts=opts.TooltipOpts(
            # trigger="item",formatter="{b}: {c} ({d}%)")
            )
        )
    return pie.dump_options_with_quotes()

@clonalType_api.route("/clonalType",methods=["POST"])
def ct():
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    if tab.shape[0] == 0:
        tab = []
        return []
    else:
        result = make_clonality(tab)
        #mycol.insert_one({"analysis":"clonalType","data":result})
        return result