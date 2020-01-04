#!/bin/bash

python tool_rpThermo.py -inputTar test_input.tar -pathway_id rp_pathway -outputTar test_output.tar

mv test_output.tar results/
