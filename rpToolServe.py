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


###################################################################################
################################## processify #####################################
###################################################################################

#hack to stop the memory leak. Indeed it seems that looping through rpFBA and the rest causes a memory leak... According to: https://github.com/opencobra/cobrapy/issues/568 there is still memory leak issues with cobrapy. looping through hundreds of models and running FBA may be the culprit

import inspect
import traceback
import signal
from functools import wraps
from multiprocessing import Process, Queue

'''
This is to deal with an error caused by segmentation fault
'''
def handler(signum, frame):
    raise OSError('Code is throwing a segmentation fault')


class Sentinel:
    pass


def processify(func):
    '''Decorator to run a function as a process.
    Be sure that every argument and the return value
    is *pickable*.
    The created process is joined, so the code does not
    run in parallel.
    '''

    def process_generator_func(q, *args, **kwargs):
        result = None
        error = None
        it = iter(func())
        while error is None and result != Sentinel:
            try:
                result = next(it)
                error = None
            except StopIteration:
                result = Sentinel
                error = None
            except Exception:
                ex_type, ex_value, tb = sys.exc_info()
                error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
                result = None
            q.put((result, error))

    def process_func(q, *args, **kwargs):
        try:
            result = func(*args, **kwargs)
        except Exception:
            ex_type, ex_value, tb = sys.exc_info()
            error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
            result = None
        else:
            error = None

        q.put((result, error))

    def wrap_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_func.__name__ = func.__name__ + 'processify_func'
        setattr(sys.modules[__name__], process_func.__name__, process_func)

        signal.signal(signal.SIGCHLD, handler) #This is to catch the segmentation error

        q = Queue()
        p = Process(target=process_func, args=[q] + list(args), kwargs=kwargs)
        p.start()
        result, error = q.get()
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

        return result

    def wrap_generator_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_generator_func.__name__ = func.__name__ + 'processify_generator_func'
        setattr(sys.modules[__name__], process_generator_func.__name__, process_generator_func)

        signal.signal(signal.SIGCHLD, handler) #This is to catch the segmentation error

        q = Queue()
        p = Process(target=process_generator_func, args=[q] + list(args), kwargs=kwargs)
        p.start()

        result = None
        error = None
        while error is None:
            result, error = q.get()
            if result == Sentinel:
                break
            yield result
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if inspect.isgeneratorfunction(func):
            return wrap_generator_func(*args, **kwargs)
        else:
            return wrap_func(*args, **kwargs)
    return wrapper


###########################################################
################## multiprocesses run #####################
###########################################################

#This is a non-deamonic multiprocessing method that can be used in combination with processify

import multiprocessing
import multiprocessing.pool


class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class nonDeamonicPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess



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
                file_name = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                rpsbml = rpSBML.rpSBML(file_name)
                rpsbml.readSBML(sbml_path)
                rpthermo.pathway_drG_prime_m(rpsbml, pathway_id)
                rpsbml.writeSBML(tmpOutputFolder)
                rpsbml = None
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    file_name += '.rpsbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))


## Less memory effecient than the _hdd method but faster
#
#
@processify
def singleThermo_multi(file_name, sbml_path, pathway_id, tmpOutputFolder):
    rpcache = rpToolCache.rpToolCache()
    rpthermo = rpThermo.rpThermo()
    rpthermo.kegg_dG = rpcache.kegg_dG
    rpthermo.cc_preprocess = rpcache.cc_preprocess
    rpsbml = rpSBML.rpSBML(file_name)
    rpsbml.readSBML(sbml_path)
    rpthermo.pathway_drG_prime_m(rpsbml, pathway_id)
    rpsbml.writeSBML(tmpOutputFolder)
    rpsbml = None

## Multiprocessing implementation of the thermodynamics package
#
#
def runThermo_multi(inputTar, outputTar, pathway_id='rp_pathway'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            #HERE SPECIFY THE NUMBER OF CORES
            #pool = multiprocessing.Pool(processes=10)
            pool = nonDeamonicPool(processes=10)
            results = []
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                file_name = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                results.append(pool.apply_async(rpthermo.singleThermo_multi, args=(file_name, sbml_path, pathway_id, tmpOutputFolder,)))
            output = [p.get() for p in results]
            logging.info(output)
            pool.close()
            pool.join()
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
def main(inputTar, outputTar, pathway_id='rp_pathway'):
    with open(inputTar, 'rb') as inputTar_bytes:
        outputTar_bytes = io.BytesIO()
        #runThermo_hdd(rpthermo, inputTar_bytes, outputTar_bytes, pathway_id)
        runThermo_multi(rpthermo, inputTar_bytes, outputTar_bytes, pathway_id)
        ########## IMPORTANT #####
        outputTar_bytes.seek(0)
        ##########################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)
