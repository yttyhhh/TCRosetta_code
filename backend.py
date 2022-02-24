from flask import Flask
from flask import request
#from flask import jsonify
from flask_cors import CORS
import pandas as pd
import json
import tempfile
#import sys
#from functools import reduce
from bson.json_util import dumps
#import numpy as np
#from operator import itemgetter
from network import network_api
from embedding import embedding_api
from public import public_api
from normal import normal_api
from diversity import diversity_api
from clonality import clonality_api
from VJUsage import VJUsage_api
from V import Vdis_api
from J import Jdis_api
from length import len_api
# from normal import normal_data_processing
from enrichment import enrichment_api
from inputParser import parser_data
from document import document_api
from clonalType import clonalType_api
from Onlysearch import onlysearch_api
from Example import example_api
# from update_test import update_test_api

app = Flask(__name__)
CORS(app)


app.register_blueprint(network_api)
app.register_blueprint(embedding_api)
app.register_blueprint(normal_api)
app.register_blueprint(public_api)
app.register_blueprint(diversity_api)
app.register_blueprint(clonality_api)
app.register_blueprint(VJUsage_api)
app.register_blueprint(Vdis_api)
app.register_blueprint(Jdis_api)
app.register_blueprint(len_api)
app.register_blueprint(enrichment_api)
app.register_blueprint(document_api)
app.register_blueprint(clonalType_api)
app.register_blueprint(onlysearch_api)
app.register_blueprint(example_api)
# app.register_blueprint(update_test_api)

@app.route('/test_api',methods=['GET'])
def hello():
    data = {'name':'wjp','age':'22'}
    return json.dumps({
        "data":data
    })

@app.route('/upload_seq', methods=["POST"])
def upload_sequence():
    req = request.json
    ss = req["sequence"]
    acid = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    sequence = []
    msg = ""
    status = 0
    st = 200
    rel_seq = ["he"]
    try:
        file_path = tempfile.NamedTemporaryFile(prefix="simpleSearch", dir="/tmp").name
        print (file_path)
        for a in ss:
            if a == "\n":
                pass
            elif a.upper() not in acid:
                msg = "Must be amino acid characters, please double check and re-enter. The separator must be one of '\\n'."
                status = 1
                st = 400
        if status == 0:
            rel_seq = []
            if "\n" in ss:
                sequence = ss.split("\n")
                for seq in sequence:
                    if len(seq)> 20 or len(seq)<8:
                        pass
                    else:
                        rel_seq.append(seq)
                msg = "Successful upload CDR3 sequence list"
                st = 200
            if "\n" not in ss:
                rel_seq.append(ss.split(" ")[0]) 
        if len(rel_seq) == 0:
            st=400
            msg="Sequence with abnormal sequence length"
        else:            
            sequence = [ i.upper() for i in rel_seq if i != ""]
            print (sequence)
            res = pd.DataFrame()
            res["AASeq"] = sequence
            #res = res.dropna()
            res.to_csv(file_path)

        return dumps(
            {
                "status":st,
                "msg":msg,
                "fipa":file_path
            }
        )
    except BaseException as e:
        print(e)
        return dumps(
            {
                "msg": "Something is wrong",
                "status": 400,
            }
        )

    
@app.route('/upload', methods=["POST"])
def handle_upload():

    try:

        file_path = tempfile.NamedTemporaryFile(prefix="TCRosseta", dir="/tmp").name

        #/tmp/TCRossetawr9tuip3，这里上传文件保存路径
        print (file_path)

        with open(file_path, 'w') as handle:
            handle.write( request.files['img'].read().decode() )

        if 'mixcr' in request.files['img'].filename.lower():
            print ("mixcr")
            res = parser_data(file_path, 'MiXCR')
        if 'catt' in request.files['img'].filename.lower():
            print ("catt")
            res = parser_data(file_path, 'CATT')
        if 'imseq' in request.files['img'].filename.lower():
            print ("imseq")
            res = parser_data(file_path, 'IMSEQ')
        if 'rtcr' in request.files['img'].filename.lower():
            print ("rtcr")
            res = parser_data(file_path, 'RTCR')
        if 'tcrosetta' in request.files['img'].filename.lower():
            print ("TCRosetta")
            res = parser_data(file_path, 'TCRosetta')
        if 'cdr3' in request.files['img'].filename.lower():
            print ("CDR3 list")
            res = parser_data(file_path,"CDR3")

        if res.shape[0] < 1000 and 'cdr3' not in request.files['img'].filename.lower():
            print (res.head())
            return dumps(
                {
                    "msg": "Low quality file.",
                    "status": 300,
                    }
                    )
        else:
        #res = pd.read_table(file_path)
            res.to_csv(file_path)
            print (res.head())
            return dumps(
                {
                    "msg": "IM OK",
                    "status": 200,
                    "filePath": file_path
                }
            )      

    except BaseException as e:

        print(e)

        return dumps(
            {
                "msg": "File Format for CATT, MIXCR, IMSEQ, is not right",
                "status": 400,
            }
        )

@app.route('/upload_example', methods=["POST"])
def example_upload():
    try:
        #/tmp/Exampleawr9tuip3，这里上传文件保存路径
        req = request.json
        filename = req["ef"]["url"]
        name = req["ef"]["name"]
        if "catt" in filename.lower():
            print ("catt")
            res = parser_data(filename, 'CATT')
        if "mixcr" in filename.lower():
            print ("mixcr")
            res = parser_data(filename, 'MiXCR')
        if "imseq" in filename.lower():
            print ("imseq")
            res = parser_data(filename, 'IMSEQ')
        if "rtcr" in filename.lower():
            print ("rtcr")
            res = parser_data(filename, 'RTCR')
        if 'tcrosetta' in filename.lower():
            print ("TCRosetta")
            res = parser_data(filename, 'TCRosetta')
        if 'cdr3' in filename.lower():
            print ("CDR3 list")
            res = parser_data(filename,"CDR3")      

        print (res.head())
        if "Vregion" in list(res.columns):
            file_path = tempfile.NamedTemporaryFile(prefix="Example", dir="/tmp").name
        else:
            file_path = tempfile.NamedTemporaryFile(prefix="ExampleEasy", dir="/tmp").name
        result = {"name":name,"url":file_path}
        res.to_csv(file_path)
        print (result)

        return dumps(
            {
                "msg": "IM OK",
                "status": 200,
                "filePath": result
            }
        )      

    except BaseException as e:

        print(e)

        return dumps(
            {
                "msg": "Something is wrong",
                "status": 400,
            }
        )

if __name__ == '__main__':
    # app.debug = True
    app.run()
