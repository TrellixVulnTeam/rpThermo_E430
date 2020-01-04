FROM brsynth/rpcache-rest

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

COPY component_contribution /home/component_contribution/
COPY input_cache.tar.xz /home/

RUN tar xf /home/input_cache.tar.xz -C /home/ && \
    rm /home/input_cache.tar.xz && \
    tar xf /home/component_contribution/component_contribution_data.tar.xz -C /home/component_contribution/ && \
    mkdir /home/cache/ && \
    mv /home/input_cache/cc_preprocess.npz /home/cache/cc_preprocess.npz
