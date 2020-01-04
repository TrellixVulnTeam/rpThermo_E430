FROM brsynth/rpcache-rest

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

COPY component_contribution /home/component_contribution/
COPY cc_preprocess.npz /home/cache/
