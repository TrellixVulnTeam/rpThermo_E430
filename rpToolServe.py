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
import shutil

sys.path.insert(0, '/home/')
import rpTool as rpThermo
import rpToolCache
import rpSBML


##
#
#
def singleThermo_mem(rpthermo, member_name, rpsbml_string, pathway_id):
    #open one of the rp SBML files
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    rpthermo.pathway_drG_prime_m(rpsbml, pathway_id)
    return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')


##
#
#
def runThermo_mem(rpthermo, inputTar, outTar, pathway_id):
    #loop through all of them and run FBA on them
    with tarfile.open(fileobj=outTar, mode='w:xz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = singleThermo_mem(rpthermo,
                            member.name,
                            in_tf.extractfile(member).read().decode("utf-8"),
                            pathway_id)
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
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
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    fileName += '.rpsbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))

##
#
#
def main(inputTar, outputTar, pathway_id='rp_pathway'):
    with open(inputTar, 'rb') as inputTar_bytes:
        outputTar_bytes = io.BytesIO()
        rpcache = rpToolCache.rpToolCache()
        rpthermo = rpThermo.rpThermo()
        rpthermo.kegg_dG = rpcache.kegg_dG
        rpthermo.cc_preprocess = rpcache.cc_preprocess
        runThermo_hdd(rpthermo, inputTar_bytes, outputTar_bytes, pathway_id)
        ########## IMPORTANT #####
        outputTar_bytes.seek(0)
        ##########################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)
