#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpThermo REST service

"""
import argparse

import rpToolServe

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rpthermo = rpThermo.rpThermo()
    rpthermo.kegg_dG = rpcache.kegg_dG
    rpthermo.cc_preprocess = rpcache.cc_preprocess
    rpToolServe.main(rpthermo, params.inputTar, params.outputTar, params.pathway_id)
