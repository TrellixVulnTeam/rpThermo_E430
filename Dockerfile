FROM brsynth/rpcache-rest:dev

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy
