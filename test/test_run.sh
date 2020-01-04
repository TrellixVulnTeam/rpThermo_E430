#!/bin/sh

docker run -d -p 8888:8888 --name test_rpThermo brsynth/rpthermo
sleep 10
python tool_rpThermo.py -inputTar test_input.tar -outputTar test_output.tar -pathway_id rp_pathway -server_url http://0.0.0.0:8888/REST 
docker kill test_rpThermo
docker rm test_rpThermo
