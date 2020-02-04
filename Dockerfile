FROM brsynth/rpcache-rest:dev

RUN conda install -y -c openbabel openbabel && \
    conda install -y -c anaconda pandas && \
    conda install -c anaconda scipy

RUN echo 'deb https://melchior.dulac@inra.fr:AKCp5ejxp21AdB9rVLTJkUrq4SoZZXNkBLgj9XTAQj7jcKV1EjrBeiV1wCce9eqvtoG6cNbPe@hub.chemaxon.com/artifactory/cxn-deb-release-local all main' >> /etc/apt/sources.list
RUN apt-get install -y apt-transport-https gnupg apt-utils ca-certificates libuuid1 libblkid-dev openjdk-8-jdk
RUN apt-get update
RUN apt-get install -y --allow-unauthenticated marvin

COPY license.cxl /home/
ENV CHEMAXON_LICENSE_URL /home/license.cxl
