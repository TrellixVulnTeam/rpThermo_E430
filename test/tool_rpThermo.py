#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpThermo REST service
python tool_rpThermo.py -inputTar test_input.tar -pathway_id rp_pathway -outputTar test_output.tar
python tool_rpThermo.py -sbml test_rpFBA.rpsbml.xml -pathway_id rp_pathway -outputTar test_output.tar


"""
import argparse
import rpToolServe
import tempfile
import tarfile

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-sbml', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    if params.sbml=='None' or params.sbml==None or params.sbml=='':
        if params.inputTar=='None' or params.inputTar==None or params.inputTar=='':	
            logging.error('Cannot have no SBML and no TAR input')
            exit(0) 
        rpToolServe.main(params.inputTar, params.outputTar, params.pathway_id)
    else:
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inT = tmpOutputFolder+'/tmp_input.tar.xz'
            with tarfile.open(inT, mode='w:xz') as tf:
                tf.add(params.sbml)
            rpToolServe.main(inT, params.outputTar, params.pathway_id)
