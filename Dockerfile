FROM biocontainers/biocontainers:debian-stretch-backports
MAINTAINER Yasset Perez-Riverol <ypriverol@gmail.com>
LABEL software="pypgatk" \
    container="pypgatk" \
    software.version="0.0.5" \
    version="1"

RUN apt-get update && apt-get install -y python3 && apt-get install -y python3-pip && apt-get install -y git
RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

WORKDIR /data
RUN mkdir -p /tool/source && cd /tool/source

RUN git config --global http.sslVerify false
RUN git clone --depth 1 https://github.com/bigbio/py-pgatk.git /tool/source
RUN cd /tool/source/ && pip3 install -r requirements.txt && python3 setup.py install

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH=$PATH:/tool/source/pypgatk/
RUN chmod +x /tool/source/pypgatk/pypgatk_cli.py
