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
import logging

sys.path.insert(0, '/home/')
import rpTool as rpThermo
import rpToolCache
import rpSBML


logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

logging.disable(logging.INFO)
logging.disable(logging.WARNING)


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
    with tarfile.open(fileobj=outTar, mode='w:gz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r') as in_tf:
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




###################### Multi ###############


## Seperate an array into equal lengths
#
#
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

### concurent

import concurrent.futures

## Multiprocessing implementation of the thermodynamics package
#
#
def runThermo_multi_concurrent(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
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
                        logging.error('%r generated an exception: %s' % (f_n, exc))
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpThermo has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    file_name += '.sbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True

### multiprocessing

import multiprocessing

## Multiprocessing implementation of the thermodynamics package
#
#
def runThermo_multi_process(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            ### construct the processes list and start
            processes = []
            for s_l in chunkIt(glob.glob(tmpInputFolder+'/*'), num_workers):
                p = multiprocessing.Process(target=singleThermo, args=(s_l, pathway_id, tmpOutputFolder,))
                processes.append(p)
                p.start()
            #wait for all to finish
            for process in processes:
                process.join()
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpThermo has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    file_name += '.sbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True

############################# single core ##########################

def runThermo_hdd(inputTar, outputTar, pathway_id='rp_pathway'):
    rpcache = rpToolCache.rpToolCache()
    rpthermo = rpThermo.rpThermo()
    rpthermo.kegg_dG = rpcache.kegg_dG
    rpthermo.cc_preprocess = rpcache.cc_preprocess
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                rpthermo.pathway_drG_prime_m(rpsbml, pathway_id)
                rpsbml.writeSBML(tmpOutputFolder)
                rpsbml = None
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpThermo has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    fileName += '.sbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True


##
#
#
def main(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway'):
    if num_workers<=0:
        logging.error('Cannot have less or 0 workers')
        return False
    elif num_workers>20:
        logging.error('20 or more is a little too many number of workers')
    elif num_workers==1:
        runThermo_hdd(inputTar, outputTar, pathway_id)
    else:
        runThermo_multi_concurrent(inputTar, outputTar, num_workers, pathway_id)
        runThermo_multi_process(inputTar, outputTar, num_workers, pathway_id)


