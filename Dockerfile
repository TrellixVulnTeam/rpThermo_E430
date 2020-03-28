FROM brsynth/rpcache:dev

RUN apt-get update
RUN apt-get install -y openjdk-8-jdk
#WARNING: the deb version of marvin must be added to the root folder
COPY marvin_linux_20.9.deb /home/
RUN dpkg -i /home/marvin_linux_20.9.deb
RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

#WARNING: make sure that you download the Marvin licence and paste it at the root of the Dockerfile
COPY license.cxl /home/
ENV CHEMAXON_LICENSE_URL /home/license.cxl

COPY component_contribution /home/component_contribution/
COPY input_cache.tar.xz /home/
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY rpToolCache.py /home/
COPY tool_rpThermo.py /home/

RUN tar xf /home/input_cache.tar.xz -C / && \
    mv /input_cache/* /home/input_cache/ && \
    rm /home/input_cache.tar.xz && \
    tar xf /home/component_contribution/component_contribution_data.tar.xz -C /home/component_contribution/ && \
    mv /home/input_cache/cc_preprocess.npz /home/cache/cc_preprocess.npz

RUN python rpToolCache.py
