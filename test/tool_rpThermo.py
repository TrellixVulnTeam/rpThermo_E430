#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpThermo REST service

"""
import requests
import argparse
import json
import tempfile
import tarfile

##
#
#
def rpThermoUpload(inputTar,
        pathway_id,
        server_url,
        outputTar):
    # Post request
    data = {'pathway_id': pathway_id}
    files = {'inputTar': open(inputTar, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server_url+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-sbml', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-server_url', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    if params.sbml=='None' or params.sbml==None or params.sbml=='':
        if params.inputTar=='None' or params.inputTar==None or params.inputTar=='':	
            logging.error('Cannot have no SBML and no TAR input')
            exit(0) 
        rpThermoUpload(params.inputTar,
                params.pathway_id,
                params.server_url,
                params.outputTar)
    else:
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inT = tmpOutputFolder+'/tmp_input.tar.xz'
            with tarfile.open(inT, mode='w:xz') as tf:
                tf.add(params.sbml)
            rpThermoUpload(inT,
                    params.pathway_id,
                    params.server_url,
                    params.outputTar)
