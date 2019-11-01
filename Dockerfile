FROM brsynth/rpbase

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy
    #conda install -y -c conda-forge openbabel && \

RUN apt-get install --quiet --yes libxext6 libxrender-dev

RUN mkdir input_cache && \
    mkdir component_contribution

COPY rpThermo.py /home/
COPY component_contribution/ component_contribution/
COPY input_cache.tar.xz /home/

#copy and extract the input_cache
#add the component contribution files
RUN tar xf /home/input_cache.tar.xz -C /home/ && \
    rm /home/input_cache.tar.xz && \
    tar xf /home/component_contribution/component_contribution_data.tar.xz -C /home/component_contribution/ && \
    mkdir /home/cache/ && \
    mv /home/input_cache/cc_preprocess.npz /home/cache/cc_preprocess.npz

RUN python /home/rpThermo.py
