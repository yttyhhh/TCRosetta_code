import math
import numpy as np
import json
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from pyecharts.globals import CurrentConfig
from pyecharts import options as opts
from pyecharts.globals import ThemeType
from pyecharts.charts import Line
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import normal_data_processing
#from Example import mycol

CurrentConfig.GLOBAL_ENV = Environment(loader=FileSystemLoader("./templates"))

diversity_api = Blueprint('diversity', __name__)

#多样性
#aplha选择不同，会导致会有不一样的情况出现，这里要设定alpha的值
#可以参照：https://github.com/vegandevs/vegan/blob/master/R/renyi.R
#编辑一个python的版本
#ll格式这里设定一下
#ll是Ratio
#renyi中的p指的是unique的TCR序对应的cloneFraction,N代表着多有unique的序列的数目
def renyi_entropy(ll,alpha = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64]):
    renyi = []
    #total 所有unique sequence的数目
    total = ll.shape[0]
    #在计算之前先要标准化处理一下
    sum_ratio = sum(ll["Ratio"])
    ll["Ratio"] = ll["Ratio"]/sum_ratio
    for a in alpha:
        if a == 0:
            renyi.append({"renyi":np.log(total),"alpha":0})
        elif a == 1:
            renyi.append({"renyi":-sum([ number * np.log(number) for number in ll["Ratio"]]),"alpha":1})
        elif a == 2:
            renyi.append({"renyi":-np.log(sum([ math.pow(number,2) for number in ll["Ratio"]])),"alpha":2})
        else:
            renyi.append({"renyi":(np.log(sum([ math.pow(number,a) for number in ll["Ratio"]])))/(1-a),"alpha":a})
    re = pd.DataFrame(renyi)
    return re

def diversity(tab):
    re = renyi_entropy(tab)
    line_x_data = list(re.alpha)
    line_y_data = list(re.renyi)
    line_x_data = [str(x) for x in line_x_data]
    line = (
        Line()
        .add_xaxis(
            xaxis_data=line_x_data,
        )
        .add_yaxis(
            series_name="Renyi index",
            y_axis=line_y_data,
            markpoint_opts=opts.MarkPointOpts(
                data=[opts.MarkPointItem(name="Shannopy entropy", coord=[line_x_data[3], line_y_data[3]])]
            ),
            is_smooth=True,
            label_opts=opts.LabelOpts(is_show=False),
            itemstyle_opts=opts.ItemStyleOpts(color="#54123B")
        )
        .set_global_opts(
            legend_opts=opts.LegendOpts(pos_top="4%"),
            xaxis_opts=opts.AxisOpts(
                type_="category",
                name = "Alpha",
                name_location="center",
                name_gap=35,
                name_textstyle_opts=opts.TextStyleOpts(font_size=20)
            ),
            yaxis_opts=opts.AxisOpts(
                type_="value",
                axistick_opts=opts.AxisTickOpts(is_show=True),
                splitline_opts=opts.SplitLineOpts(is_show=True),
                name = "Renyi Index",
                name_location="center",
                name_gap=35,
                name_textstyle_opts=opts.TextStyleOpts(font_size=20)
            ),
            toolbox_opts = opts.ToolboxOpts(is_show=True,
                                    pos_top="3%",
                                    pos_right="right",
                                    feature={
                                        "saveAsImage":{
                                            "title":"save as image"
                                        }
                                    })
        )
    )
    return line.dump_options_with_quotes()

@diversity_api.route("/diversity",methods=["POST","GET"])
def vj():
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    if tab.shape[0] == 0:
        tab = []
        return tab
    else:
        result = diversity(tab)
        #mycol.insert_one({"analysis":"diveristy","data":result})
        return result