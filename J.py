import pandas as pd
from jinja2 import Environment, FileSystemLoader
from pyecharts.globals import CurrentConfig
from pyecharts import options as opts
from pyecharts.charts import Bar
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import normal_data_processing
from Example import mycol

CurrentConfig.GLOBAL_ENV = Environment(loader=FileSystemLoader("./templates"))

Jdis_api = Blueprint('J', __name__)

j_dis_background = pd.read_csv("/workspace2/yuet/tcrosetta/background_data/background_dis_j.csv")

def j_dis(tab):
    tab["Count_j"] = tab.groupby("Jregion")["Condition"].transform("count")
    tab = tab.sort_values(by=["Count_j"],ascending=False)
    tab_j = tab.drop_duplicates(subset=['Jregion'],keep='first')
    j_dis = tab_j[["Jregion","Count_j"]]
    j_dis_result = pd.merge(j_dis,j_dis_background)
    j_dis_result.fillna(0,inplace=True)
    x = [str(i) for i in j_dis_result["Jregion"]]
    y1 = [int(i) for i in j_dis_result["Count_j"]]
    y2 = [int(i) for i in j_dis_result["Count_j_bg"]]
    bar = (
        Bar()
        .add_xaxis(x)
        .add_yaxis("Sample",y1,itemstyle_opts=opts.ItemStyleOpts(color="#54123B"))
        .add_yaxis("Background",y2, is_selected=False,itemstyle_opts=opts.ItemStyleOpts(color="#81382d"))
        .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
        .set_global_opts(
            legend_opts=opts.LegendOpts(pos_top="4%"),
            title_opts=opts.TitleOpts(title=""),
            xaxis_opts=opts.AxisOpts(
                axislabel_opts=opts.LabelOpts(rotate=70)),
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

@Jdis_api.route("/Jdis",methods=["POST"])
def jd():
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    if tab.shape[0] == 0:
        tab = []
        return tab
    else:
        result = j_dis(tab)
        #mycol.insert_one({"analysis":"J","data":result})
        return result