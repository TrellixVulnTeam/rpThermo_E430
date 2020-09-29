#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Extract the sink from an SBML into RP2 friendly format

"""
import argparse
import tempfile
import os
import logging
import shutil
import docker


##
#
#
def main(inputfile, 
         input_format,
         output,
         pathway_id='rp_pathway',
         ph=7.0,
         ionic_strength=200.0,
         pMg=10.0,
         temp_k=298.15):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpthermo-standalone:equilibrator'
    try:
        image = docker_client.images.get(image_str)
    except docker.errors.ImageNotFound:
        logging.warning('Could not find the image, trying to pull it')
        try:
            docker_client.images.pull(image_str)
            image = docker_client.images.get(image_str)
        except docker.errors.ImageNotFound:
            logging.error('Cannot pull image: '+str(image_str))
            exit(1)
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        shutil.copy(inputfile, tmpOutputFolder+'/input.dat')
        command = ['/home/tool_rpThermo.py',
                   '-input',
                   '/home/tmp_output/input.dat',
                   '-input_format',
                   str(input_format),
                   '-pathway_id',
                   str(pathway_id),
                   '-ph',
                   str(ph),
                   '-ionic_strength',
                   str(ionic_strength),
                   '-pMg',
                   str(pMg),
                   '-temp_k',
                   str(temp_k),
                   '-output',
                   '/home/tmp_output/output.dat']
        container = docker_client.containers.run(image_str, 
                                                 command, 
                                                 detach=True, 
                                                 stderr=True, 
                                                 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=False, stderr=True)
        err_str = err.decode('utf-8')
        if not 'ERROR' in err_str:
            shutil.copy(tmpOutputFolder+'/output.dat', output)
        else:
            print(err_str)
        container.remove()




##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-input', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-ph', type=float, default=7.0)
    parser.add_argument('-ionic_strength', type=float, default=200.0)
    parser.add_argument('-pMg', type=float, default=10.0)
    parser.add_argument('-temp_k', type=float, default=298.15)
    params = parser.parse_args()
    main(params.input, params.input_format, params.output, params.pathway_id, params.ph, params.ionic_strength, params.pMg, params.temp_k)
