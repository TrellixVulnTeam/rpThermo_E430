#FROM brsynth/rpcache-rest
FROM brsynth/rpcache

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

COPY component_contribution /home/component_contribution/
COPY input_cache.tar.xz /home/
COPY license.cxl /home/
ENV CHEMAXON_LICENSE_URL /home/license.cxl
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY rpToolCache.py /home/

RUN tar xf /home/input_cache.tar.xz -C / && \
    mv /input_cache/* /home/input_cache/ && \
    rm /home/input_cache.tar.xz && \
    tar xf /home/component_contribution/component_contribution_data.tar.xz -C /home/component_contribution/ && \
    #mkdir /home/cache/ && \
    mv /home/input_cache/cc_preprocess.npz /home/cache/cc_preprocess.npz

RUN python rpToolCache.py
