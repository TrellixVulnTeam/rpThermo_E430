#!/bin/bash

docker run -v ${PWD}/inside_run_sbml.sh:/home/inside_run_sbml.sh -v ${PWD}/tool_rpThermo.py:/home/tool_rpThermo.py -v ${PWD}/test_rpCofactors.rpsbml.xml:/home/test_rpCofactors.rpsbml.xml -v ${PWD}/results/:/home/results/ --rm brsynth/rpthermo /bin/sh /home/inside_run_sbml.sh

cp results/test_output.tar .
