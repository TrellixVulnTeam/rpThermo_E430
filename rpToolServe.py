#!/usr/bin/env python3

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


import concurrent.futures

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



def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out


## Less memory effecient than the _hdd method but faster
#
#
def singleThermo(sbml_paths, pathway_id, tmpOutputFolder):
    rpcache = rpToolCache.rpToolCache()
    rpthermo = rpThermo.rpThermo()
    rpthermo.kegg_dG = rpcache.kegg_dG
    rpthermo.cc_preprocess = rpcache.cc_preprocess
    for sbml_path in sbml_paths:
        file_name = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
        rpsbml = rpSBML.rpSBML(file_name)
        rpsbml.readSBML(sbml_path)
        rpthermo.pathway_drG_prime_m(rpsbml, pathway_id)
        rpsbml.writeSBML(tmpOutputFolder)
        rpsbml = None


## Multiprocessing implementation of the thermodynamics package
#
#
def runThermo_multi(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
                jobs = {}
                #split the files "equally" between all workers
                for s_l in chunkIt(glob.glob(tmpInputFolder+'/*'), num_workers):
                    jobs[executor.submit(singleThermo, s_l, pathway_id, tmpOutputFolder)] = s_l
                for future in concurrent.futures.as_completed(jobs):
                    f_n = jobs[future]
                    try:
                        data = future.result()
                    except Exception as exc:
                        print('%r generated an exception: %s' % (f_n, exc))
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    file_name += '.rpsbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))



##
#
#
def main(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway'):
    with open(inputTar, 'rb') as inputTar_bytes:
        outputTar_bytes = io.BytesIO()
        runThermo_multi(inputTar_bytes, outputTar_bytes, num_workers, pathway_id)
        ########## IMPORTANT #####
        outputTar_bytes.seek(0)
        ##########################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)
