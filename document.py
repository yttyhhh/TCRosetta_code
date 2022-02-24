import pymongo
from flask_cors import CORS
from flask import Flask,request, Blueprint

document_api = Blueprint('document', __name__)

#这里1和10的mongo都再用
myclient = pymongo.MongoClient(
    "mongodb://222.20.95.98:26906",
    username="yuet",
    password="HUST11405yt323.",
    authMechanism="SCRAM-SHA-1")
mydb = myclient["document"]
mycol = mydb["static"]
mycol3 = mydb["static2"]


myclient2 = pymongo.MongoClient(
    "mongodb://222.20.95.102:26906",
    username="yuet",
    password="HUST11405yt323.",
    authMechanism="SCRAM-SHA-1")
# myclient2 = pymongo.MongoClient("mongodb://222.20.95.102:26906")
mydb2 = myclient2["document"]
mycol2 = mydb2["home_static"]


print (1)

#document后面document可能变成图片的形式
@document_api.route("/document_length",methods=["GET"])
def le():
    length = list(mycol3.find({"cond":"length"},{"data"}))[0]["data"]
    return length
@document_api.route("/document_vj",methods=["GET"])
def v_j():
    vj = list(mycol3.find({"cond":"vj"},{"data"}))[0]["data"]
    return vj
@document_api.route("/document_diversity",methods=["GET"])
def div():
    diversity = list(mycol3.find({"cond":"diversity"},{"data"}))[0]["data"]
    return diversity
@document_api.route("/document_clonality",methods=["GET"])
def clo():
    clonality = list(mycol3.find({"cond":"clonality"},{"data"}))[0]["data"]
    return clonality
@document_api.route("/document_public",methods=["GET"])
def pu():
    public = list(mycol.find({"cond":"public"},{"data"}))[0]["data"]
    return public
@document_api.route("/document_network",methods=["GET"])
def net():
    network = list(mycol.find({"cond":"network"},{"data"}))[0]["data"]
    return network
@document_api.route("/document_enrichment",methods=["GET"])
def en():
    enrichment = list(mycol.find({"cond":"enrichment"},{"data"}))[0]["data"]
    return enrichment
@document_api.route("/document_v",methods=["GET"])
def v_dis():
    v = list(mycol3.find({"cond":"v"},{"data"}))[0]["data"]
    return v
@document_api.route("/document_j",methods=["GET"])
def j_dis():
    j = list(mycol3.find({"cond":"j"},{"data"}))[0]["data"]
    return j
@document_api.route("/document_clonaltype",methods=["GET"])
def ctype():
    clt = list(mycol3.find({"cond":"clonaltype"},{"data"}))[0]["data"]
    return clt


#这里是home页面的走马灯
@document_api.route("/home_length",methods=["GET"])
def hle():
    length = list(mycol2.find({"cond":"length"},{"data"}))[0]["data"]
    return length
@document_api.route("/home_vj",methods=["GET"])
def hv_j():
    vj = list(mycol2.find({"cond":"vj"},{"data"}))[0]["data"]
    return vj
@document_api.route("/home_network",methods=["GET"])
def hnet():
    network = list(mycol2.find({"cond":"network"},{"data"}))[0]["data"]
    return network
@document_api.route("/home_clonaltype",methods=["GET"])
def hctype():
    clt = list(mycol2.find({"cond":"clonaltype"},{"data"}))[0]["data"]
    return clt
@document_api.route("/home_clonality",methods=["GET"])
def hclo():
    clonality = list(mycol2.find({"cond":"clonality"},{"data"}))[0]["data"]
    return clonality

