#FROM brsynth/rpcache-rest
FROM brsynth/rpcache

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

COPY rpTool.py /home/
COPY rpToolServe.py /home/
