import pandas as pd
import numpy as np
from jinja2 import Environment, FileSystemLoader
from pyecharts.globals import CurrentConfig
from pyecharts import options as opts
from pyecharts.globals import ThemeType
from pyecharts.charts import Bar3D
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import normal_data_processing
from Example import mycol

CurrentConfig.GLOBAL_ENV = Environment(loader=FileSystemLoader("./templates"))

VJUsage_api = Blueprint('VJUsage', __name__)

#VJl联合使用
def vjusage(tab):
    tab = tab[["Vregion","Jregion","Ratio"]]
    vj = tab.pivot_table(index="Vregion",columns="Jregion",values="Ratio",aggfunc=np.sum)
    vj.fillna(0,inplace=True)
    j_col = list(vj.columns.values)
    v_ind = list(vj.index.values)
    list2 = []
    for idx,ins in enumerate(v_ind):
        for m,n in enumerate(j_col):
            list1 = [ins,n,vj.at[ins,n]]
            list2.append(list1)
    max_data = max([ d[2] for d in list2])
    bar3d = (
        Bar3D(init_opts=opts.InitOpts())
            .add(
                series_name="",
                data=list2,
                xaxis3d_opts=opts.Axis3DOpts(
                    type_="category", 
                    data=v_ind,
                    name="V gene",
                    name_gap = 23,
                    textstyle_opts=opts.TextStyleOpts(font_size=17)),
                yaxis3d_opts=opts.Axis3DOpts(
                    type_="category", 
                    data=j_col,
                    name="J gene",
                    name_gap = 23,
                    textstyle_opts=opts.TextStyleOpts(font_size=17)),
                zaxis3d_opts=opts.Axis3DOpts(
                    type_="value",
                    name="Frequency",
                    name_gap = 25,
                    textstyle_opts=opts.TextStyleOpts(font_size=17)),
            )
            .set_global_opts(
                visualmap_opts=opts.VisualMapOpts(
                    max_=max_data,
                    range_color=[
                        "#313695",
                        "#4575b4",
                        "#74add1",
                        "#abd9e9",
                        "#e0f3f8",
                        "#ffffbf",
                        "#fee090",
                        "#fdae61",
                        "#A05F96",
                        "#971549",
                        "#54123B",
                    ],
                ),
                toolbox_opts=opts.ToolboxOpts(is_show=True,
                                        pos_top="3%",
                                        pos_right="right",
                                        feature={
                                            "saveAsImage":{
                                                "title":"save as image"
                                            }
                                        })
                                        )
            )
    return bar3d.dump_options_with_quotes()

@VJUsage_api.route("/vjusage",methods=["POST"])
def vj():
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    if tab.shape[0] == 0:
        tab = []
        return tab
    else:
        result = vjusage(tab)
        #mycol.insert_one({"analysis":"VJ","data":result})
        return result
