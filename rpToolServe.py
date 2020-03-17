#!/usr/bin/env python3

#from contextlib import closing
#import time
import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os

import io
#import zipfile
import tarfile
import glob
import tempfile

import json
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api

sys.path.insert(0, '/home/')
import rpTool as rpThermo
import rpToolCache
import rpSBML


##
#
#
def runThermo_mem(rpthermo, inputTar, outputTar, pathway_id):
    #loop through all of them and run FBA on them
    with tarfile.open(fileobj=outputTar, mode='w:xz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    fileName = member.name.replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                    rpsbml = rpSBML.rpSBML(filename, libsbml.readSBMLFromString(in_tf.extractfile(member).read().decode("utf-8")))
                    rpthermo.pathway_drG_prime_m(rpsbml, pathway_id)
                    sbml_bytes = libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
                    fiOut = io.BytesIO(sbml_bytes)
                    info = tarfile.TarInfo(fileName+'.rpsbml.xml')
                    info.size = len(sbml_bytes)
                    tf.addfile(tarinfo=info, fileobj=fiOut)


##
#
#
def runThermo_hdd(rpthermo, inputTar, outputTar, pathway_id='rp_pathway'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                rpthermo.pathway_drG_prime_m(rpsbml, pathway_id)
                rpsbml.writeSBML(tmpOutputFolder)
                rpsbml = None
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))+'.rpsbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))


#######################################################
############## REST ###################################
#######################################################


app = Flask(__name__)
api = Api(app)


#global thermo parameter
rpcache = rpToolCache.rpToolCache()

def stamp(data, status=1):
    appinfo = {'app': 'rpThermo', 'version': '1.0',
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(),
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


class RestApp(Resource):
    """ REST App."""
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


class RestQuery(Resource):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):
        inputTar = request.files['inputTar']
        params = json.load(request.files['data'])
        #pass the files to the rpReader
        outputTar = io.BytesIO()
        rpthermo = rpThermo.rpThermo()
        rpthermo.kegg_dG = rpcache.kegg_dG
        rpthermo.cc_preprocess = rpcache.cc_preprocess
        ###### MEM ######
        #runThermo_mem(rpthermo, inputTar, outputTar, params['pathway_id'])
        ###### HDD ######
        runThermo_hdd(rpthermo, inputTar, outputTar, params['pathway_id'])
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpThermo.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=True)
