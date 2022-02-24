import pandas as pd
from collections import Counter
from jinja2 import Environment, FileSystemLoader
from pyecharts.globals import CurrentConfig
from pyecharts import options as opts
from pyecharts.charts import Bar
from flask_cors import CORS
from flask import Flask,request, Blueprint
from pyecharts.options.global_options import AxisOpts
from pyecharts.types import Axis
from normal import normal_data_processing
from Example import mycol

CurrentConfig.GLOBAL_ENV = Environment(loader=FileSystemLoader("./templates"))

len_api = Blueprint('length', __name__)

bg_length = pd.read_csv("/workspace2/yuet/tcrosetta/background_data/background_length_dis.csv")


def len_dis(tab):
    sample_length = [ len(ins["AASeq"]) for idx,ins in tab.iterrows()]
    sl = [ {"Length":idx,"Number":ins} for idx,ins in Counter(sample_length).items()]
    sl = pd.DataFrame(sl)
    len_dis_result = pd.merge(sl,bg_length)
    len_dis_result = len_dis_result.fillna(0)
    len_dis_result.sort_values(by=["Length"],ascending=True,inplace=True)
    x = [str(i) for i in len_dis_result["Length"]]
    y1 = [int(i) for i in len_dis_result["Number"]]
    y2 = [int(i) for i in len_dis_result["Number_bg"]]
    bar = (
        Bar()
        .add_xaxis(x)
        .add_yaxis("Sample",y1,itemstyle_opts=opts.ItemStyleOpts(color="#54123B"))
        .add_yaxis("Background",y2, is_selected=False,itemstyle_opts=opts.ItemStyleOpts(color="#81382d"))
        .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
        .set_global_opts(
            legend_opts=opts.LegendOpts(pos_top="4%"),
            title_opts=opts.TitleOpts(),
            xaxis_opts=opts.AxisOpts(
                axislabel_opts=opts.LabelOpts(rotate=70),
                name = "Length",
                name_location="center",
                name_gap=30,
                name_textstyle_opts=opts.TextStyleOpts(font_size=20)
                ),
            yaxis_opts = opts.AxisOpts(
                name = "Number",
                name_location="center",
                name_gap=35,
                name_textstyle_opts=opts.TextStyleOpts(font_size=20)
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
    return bar.dump_options_with_quotes()

@len_api.route("/lendis",methods=["POST"])
def ld():
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    #这里都得加个判断，判断是否为空，如果是空，就不返回图片了，返回一段话，说起质量不高
    if tab.shape[0] == 0:
        tab = []
        return tab
    else:
        result = len_dis(tab)
        #mycol.insert_one({"analysis":"length","data":result})
        return result