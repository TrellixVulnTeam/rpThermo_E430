#!/bin/bash

python tool_rpThermo.py -sbml test_rpCofactors.rpsbml.xml -pathway_id rp_pathway -outputTar test_output.tar

mv test_output.tar results/
