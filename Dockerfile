FROM brsynth/rpcache:v2

RUN apt-get update
RUN apt-get install -y openjdk-8-jdk
#WARNING: the deb version of marvin must be added to the root folder
COPY marvin_linux_20.9.deb /home/
RUN dpkg -i /home/marvin_linux_20.9.deb
RUN conda install -y -c anaconda pandas && \
    conda install -y -c anaconda scipy && \
    conda install -y -c conda-forge openbabel

#install the develop version of equilibrator
#RUN pip install equilibrator-api equilibrator-cache equilibrator-pathway
RUN git clone --single-branch --branch develop https://gitlab.com/equilibrator/equilibrator-api.git
RUN cd equilibrator-api && pip install -e . && cd ..

#equilibrator-assets
RUN git clone https://gitlab.com/equilibrator/equilibrator-assets.git
RUN cd equilibrator-assets && pip install -e . && cd ..

#equilibrator-pathway
RUN pip install equilibrator-pathway==0.3.1
#RUN git clone --single-branch --branch develop https://gitlab.com/equilibrator/equilibrator-pathway.git
#RUN cd equilibrator-pathway && pip install -e . && cd ..


#WARNING: make sure that you download the Marvin licence and paste it at the root of the Dockerfile
COPY license.cxl /home/
ENV CHEMAXON_LICENSE_URL /home/license.cxl

COPY component_contribution_legacy /home/component_contribution_legacy/
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rpThermo.py /home/
COPY galaxy/code/tool_rpMDF.py /home/
#COPY rpComponentContribution.py /home/ #if using legacy component contribution
COPY rpEquilibrator.py /home/

RUN mkdir /home/data/
COPY data/mnx_default_conc.json /home/data/

#COPY component_contribution_data.tar.xz /home/component_contribution/
RUN tar xf /home/component_contribution_legacy/component_contribution_data.tar.xz -C /home/component_contribution_legacy/
RUN rm /home/component_contribution_legacy/component_contribution_data.tar.xz

#run to downloads the required chuncks
RUN python rpEquilibrator.py
