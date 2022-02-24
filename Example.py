import pymongo
from flask_cors import CORS
from flask import Flask,request, Blueprint
import json

example_api = Blueprint('example', __name__)

myclient = pymongo.MongoClient(
    "mongodb://222.20.95.98:26906",
    username="yuet",
    password="HUST11405yt323.",
    authMechanism="SCRAM-SHA-1")
mydb = myclient["Fast_example"]
mycol = mydb["Example_data"]


@example_api.route("/example_analysis",methods=["GET"])
def example():
    print ("editor !!!!!!")
    embedding = list(mycol.find({"analysis":"embedding"},{"data":1}))[0]["data"]
    public = list(mycol.find({"analysis":"public"},{"data":1}))[0]["data"]
    enrichment = list(mycol.find({"analysis":"enrichment"},{"data":1}))[0]["data"]
    network = list(mycol.find({"analysis":"network"},{"data":1}))[0]["data"]
    onlysearch = list(mycol.find({"analysis":"onlysearch"},{"data":1}))[0]["data"]
    result = {
        "embedding":embedding,
        "public":public,
        "network":network,
        "onlysearch":onlysearch,
        "enrichment":enrichment
    }
    return json.dumps(result)
@example_api.route("/example_length",methods=["GET"])
def le():
    length = list(mycol.find({"analysis":"length"},{"data"}))[0]["data"]
    return length
@example_api.route("/example_vj",methods=["GET"])
def v_j():
    vj = list(mycol.find({"analysis":"VJ"},{"data"}))[0]["data"]
    return vj
@example_api.route("/example_diversity",methods=["GET"])
def div():
    diversity = list(mycol.find({"analysis":"diveristy"},{"data"}))[0]["data"]
    return diversity
@example_api.route("/example_clonality",methods=["GET"])
def clo():
    clonality = list(mycol.find({"analysis":"clonality"},{"data"}))[0]["data"]
    return clonality
@example_api.route("/example_V",methods=["GET"])
def V():
    vdis = list(mycol.find({"analysis":"V"},{"data"}))[0]["data"]
    return vdis
@example_api.route("/example_J",methods=["GET"])
def J():
    jdis = list(mycol.find({"analysis":"J"},{"data"}))[0]["data"]
    return jdis
@example_api.route("/example_logo",methods=["GET"])
def logo():
    lg = list(mycol.find({"analysis":"logo"},{"data"}))[0]["data"]
    return lg
@example_api.route("/example_clonaType",methods=["GET"])
def clonaType():
    ct = list(mycol.find({"analysis":"clonalType"},{"data"}))[0]["data"]
    return ct