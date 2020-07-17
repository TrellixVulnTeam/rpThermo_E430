FROM brsynth/rpcache:v1

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
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rpThermo.py /home/

#COPY component_contribution_data.tar.xz /home/component_contribution/
RUN tar xf /home/component_contribution/component_contribution_data.tar.xz -C /home/component_contribution/
RUN rm /home/component_contribution/component_contribution_data.tar.xz
