
import math
import json
import numpy as np
import pandas as pd
from pyecharts.commons.utils import JsCode
from jinja2 import Environment, FileSystemLoader
from pyecharts.globals import CurrentConfig
from pyecharts import options as opts
from pyecharts.globals import ThemeType
from pyecharts.charts import Boxplot,Grid
from flask_cors import CORS
from flask import Flask,request, Blueprint
from normal import normal_data_processing 
import pymongo
#from Example import mycol
from Example import myclient

CurrentConfig.GLOBAL_ENV = Environment(loader=FileSystemLoader("./templates"))

# myclient = pymongo.MongoClient("mongodb://222.20.95.98:26906")
mydb_clonality = myclient["tcrosetta"]
# mycol_clonality = mydb["clonality"]

clonality_api = Blueprint('clonality', __name__)

#克隆性

def clonality(ll):
    reference_clonality = list(mydb_clonality["reference"].find({"analysis":"clonality"}))[0]["data"]
    reference_xdata = list(mydb_clonality["reference"].find({"analysis":"clonality"}))[0]["xdata"]
    sum_ratio = sum(ll["Ratio"])
    ll["Ratio"] = ll["Ratio"]/sum_ratio
    total = ll.shape[0]
    all_number = []
    for number in ll["Ratio"]/np.log(total):
        rn = round(number,5)
        if rn < 0.00001:
            pass
        else:
            all_number.append(rn*np.log10(rn))
    clo = 1-(-sum(all_number))
    #clo = 1-(-sum([ number*np.log(number) for number in ll["Ratio"]])/np.log(total))
    clo = round(clo,5)
    xdata = reference_xdata
    data = reference_clonality
    xdata.insert(0,"Sample")
    data.insert(0,[clo,clo,clo,clo])
    boxlplot = Boxplot()
    final_box = boxlplot.prepare_data(data)
    data1 = [{"axis_Data":xdata,"boxData":final_box}]
    box = (
        Boxplot(init_opts=opts.InitOpts(width="800px", height="600px"))
        .add_xaxis(xaxis_data=data1[0]["axis_Data"])
        .add_yaxis(
            series_name="",
            y_axis=data1[0]["boxData"],
            # tooltip_opts=opts.TooltipOpts(
            #     formatter=JsCode(
            #         """function(param) { return [
            #                     'Disease ' + param.name + ': ',
            #                     'upper: ' + param.data[0],
            #                     'Q1: ' + param.data[1],
            #                     'median: ' + param.data[2],
            #                     'Q3: ' + param.data[3],
            #                     'lower: ' + param.data[4]
            #                 ].join('<br/>') }"""
            #     )
            # ),
        )
        .set_global_opts(
            toolbox_opts = opts.ToolboxOpts(
                is_show=True,
                pos_top="3%",
                pos_right="right",
                feature={
                    "saveAsImage":{
                        "title":"save as image" 
                    }
                }),
            # title_opts=opts.TitleOpts(
            #     title="upper: Q3 + 1.5 * IQR \nlower: Q1 - 1.5 * IQR", 
            #     pos_left="10%",
            #     pos_top="80%",
            #     title_textstyle_opts=opts.TextStyleOpts(
            #         border_color="#999", border_width=1, font_size=14
            #     )),
            legend_opts=opts.LegendOpts(pos_top="3%"),
            tooltip_opts=opts.TooltipOpts(trigger="item", axis_pointer_type="shadow"),
            xaxis_opts=opts.AxisOpts(
                name = "Disease",
                name_location="center",
                name_gap=30,
                name_textstyle_opts=opts.TextStyleOpts(font_size=20),
                boundary_gap=True,
                splitarea_opts=opts.SplitAreaOpts(
                    areastyle_opts=opts.AreaStyleOpts(opacity=1)
                ),
                axislabel_opts=opts.LabelOpts(
                    formatter="{value}",rotate=40),
                splitline_opts=opts.SplitLineOpts(is_show=False),
            ),
            yaxis_opts=opts.AxisOpts(
                type_="value",
                splitarea_opts=opts.SplitAreaOpts(is_show=False),
                name = "Clonality",
                name_location="center",
                name_gap=35,
                name_textstyle_opts=opts.TextStyleOpts(font_size=20)
                
            ),
            datazoom_opts=[
                opts.DataZoomOpts(type_="inside", range_start=0, range_end=100),
                opts.DataZoomOpts(type_="slider", xaxis_index=0, is_show=True),
            ],
        )
    )
    grid = (
        Grid()
        .add(box,grid_opts=opts.GridOpts(pos_bottom="24%"))
    )
    return grid.dump_options_with_quotes()

@clonality_api.route("/clonality",methods=['POST'])
def cl():
    req = request.json
    fp = req["fp"]
    tab = pd.read_csv(fp)
    tab = normal_data_processing(tab)
    if tab.shape[0] == 0:
        tab = []
        return tab
    else:
        result = clonality(tab)
        #mycol.insert_one({"analysis":"clonality","data":result})
        return result

