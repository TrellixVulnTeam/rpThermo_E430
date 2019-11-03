import os
import urllib.request
import pickle
import csv
import gzip
import json
import logging
import numpy as np


## Class to generate the cache
#
# Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the
#the other steps. These should be called only when the files have changes
class rpCache:
    ## Cache constructor
    #
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    def __init__(self):
        #given by Thomas
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCache')
        self.convertMNXM = {'MNXM162231': 'MNXM6',
                'MNXM84': 'MNXM15',
                'MNXM96410': 'MNXM14',
                'MNXM114062': 'MNXM3',
                'MNXM145523': 'MNXM57',
                'MNXM57425': 'MNXM9',
                'MNXM137': 'MNXM588022'}
        self.kegg_dG = None
        self.cc_preprocess = None
        if not self._loadCache():
            raise ValueError


    ##
    #
    #
    def _loadCache(self, fetchInputFiles=False):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        ###################### Fetch the files if necessary ######################
        #cc_compounds.json.gz
        if not os.path.isfile(dirname+'/input_cache/cc_compounds.json.gz') or fetchInputFiles:
            urllib.request.urlretrieve('TODO',
                    dirname+'/input_cache/cc_compounds.json.gz')
        #alberty.json
        if not os.path.isfile(dirname+'/input_cache/alberty.json') or fetchInputFiles:
            urllib.request.urlretrieve('TODO',
                    dirname+'/input_cache/alberty.json')
        #compounds.csv
        if not os.path.isfile(dirname+'/input_cache/compounds.csv') or fetchInputFiles:
            urllib.request.urlretrieve('TODO',
                    dirname+'/input_cache/compounds.csv')
        #cc_preprocess.npz
        if not os.path.isfile(dirname+'/cache/cc_preprocess.npz') or fetchInputFiles:
            urllib.request.urlretrieve('TODO',
                    dirname+'/cache/compounds.csv')
        ###################### Populate the cache #################################
        #kegg_dG.pickle
        if not os.path.isfile(dirname+'/cache/kegg_dG.pickle'):
            pickle.dump(self.kegg_dG(dirname+'/input_cache/cc_compounds.json.gz',
                dirname+'/input_cache/alberty.json',
                dirname+'/input_cache/compounds.csv'),
                open(dirname+'/cache/kegg_dG.pickle', 'wb'))
        self.kegg_dG = pickle.load(open(dirname+'/cache/kegg_dG.pickle', 'rb'))
        self.cc_preprocess = np.load(dirname+'/cache/cc_preprocess.npz')
        return True



    ## Function exctract the dG of components
    #
    #  This function parses a file from component analysis and uses KEGG id's. This means that to use precalculated
    # values in this file a given ID must have a KEGG id listed here
    #
    #  @param self Object pointer
    #  @param cc_compounds_path cc_compounds.json.gz file path
    #  @param alberty_path alberty.json file path
    #  @param compounds_path compounds.csv file path
    def kegg_dG(self,
                cc_compounds_path,
                alberty_path,
                compounds_path):
        cc_alberty = {}
        ########################## compounds ##################
        #contains the p_kas and molecule decomposition
        cid_comp = {}
        with open(compounds_path) as f:
            c = csv.reader(f, delimiter=',', quotechar='"')
            next(c)
            for row in c:
                cid_comp[row[-1].split(':')[1]] = {}
                cid_comp[row[-1].split(':')[1]]['atom_bag'] = literal_eval(row[3])
                cid_comp[row[-1].split(':')[1]]['p_kas'] = literal_eval(row[4])
                cid_comp[row[-1].split(':')[1]]['major_ms'] = int(literal_eval(row[6]))
                cid_comp[row[-1].split(':')[1]]['number_of_protons'] = literal_eval(row[7])
                cid_comp[row[-1].split(':')[1]]['charges'] = literal_eval(row[8])
        ####################### cc_compounds ############
        #TODO: seems like the new version of equilibrator got rid of this file... need to update the function
        #to take as input the new file --> i.e. the JSON input
        #notFound_cc = []
        gz_file = gzip.open(cc_compounds_path, 'rb')
        f_c = gz_file.read()
        c = json.loads(f_c)
        for cd in c:
            #find the compound descriptions
            try:
                cd.update(cid_comp[cd['CID']])
            except KeyError:
                pass
            #add the CID
            #if not mnxm in cc_alberty:
            if not cd['CID'] in cc_alberty:
                cc_alberty[cd['CID']] = {}
                if not 'alberty' in cc_alberty[cd['CID']]:
                    cc_alberty[cd['CID']]['alberty'] = [cd]
                else:
                    cc_alberty[cd['CID']]['alberty'].append(cd)
        return cc_alberty


if __name__== "__main__":
    rpcache = rpCache()
