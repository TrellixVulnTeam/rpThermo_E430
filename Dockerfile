FROM brsynth/rpcache-rest

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

RUN apt-get install --quiet --yes libxext6 libxrender-dev

COPY component_contribution/ /home/component_contribution/
