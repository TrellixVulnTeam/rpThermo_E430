#FROM brsynth/rpcache-rest
FROM brsynth/rpcache

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

RUN echo 'deb https://melchior.dulac@inra.fr:AKCp5ejxp21AdB9rVLTJkUrq4SoZZXNkBLgj9XTAQj7jcKV1EjrBeiV1wCce9eqvtoG6cNbPe@hub.chemaxon.com/artifactory/cxn-deb-release-local all main' >> /etc/apt/sources.list
RUN apt-get install -y apt-transport-https gnupg apt-utils ca-certificates libuuid1 libblkid-dev openjdk-8-jdk
RUN apt-get update
RUN apt-get install -y --allow-unauthenticated marvin

COPY component_contribution /home/component_contribution/
COPY input_cache.tar.xz /home/
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY rpToolCache.py /home/
COPY license.cxl /home/
ENV CHEMAXON_LICENSE_URL /home/license.cxl

RUN tar xf /home/input_cache.tar.xz -C / && \
    mv /input_cache/* /home/input_cache/ && \
    rm /home/input_cache.tar.xz && \
    tar xf /home/component_contribution/component_contribution_data.tar.xz -C /home/component_contribution/ && \
    #mkdir /home/cache/ && \
    mv /home/input_cache/cc_preprocess.npz /home/cache/cc_preprocess.npz

RUN python rpToolCache.py
