FROM brsynth/rpcache


RUN apt-get update
RUN apt-get install -y apt-transport-https gnupg apt-utils ca-certificates libuuid1 libblkid-dev libblkid1 openjdk-8-jdk
#RUN echo 'deb https://<E-MAIL>:<API-KEY>@hub.chemaxon.com/artifactory/cxn-deb-release-local all main' >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y --allow-unauthenticated marvin
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
