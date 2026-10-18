# syntax=docker/dockerfile:1
FROM docker.io/continuumio/miniconda3:25.3.1-1

LABEL description="HLA-HD typing container (bowtie2 + HLA-HD)" \
      author="anode"

USER root
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    wget \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
    

ARG ENV_NAME="hla_typing"

RUN conda create -y -n ${ENV_NAME} python=3.11 \
    && conda install -y -n ${ENV_NAME} -c bioconda -c conda-forge bowtie2 \
    && conda clean -afy


ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:/opt/conda/bin:$PATH

WORKDIR /opt

COPY hlahd.1.7.0.tar.gz /opt/

RUN tar -xzf /opt/hlahd.1.7.0.tar.gz -C /opt/ \
    && rm /opt/hlahd.1.7.0.tar.gz \
    && cd /opt/hlahd.1.7.0 \
    && chmod +x install.sh update.dictionary.sh \
    && bash install.sh

ENV PATH=/opt/hlahd.1.7.0/bin:$PATH
ENV HLAHD_BASE=/opt/hlahd.1.7.0

COPY run_hlahd.sh /usr/local/bin/run_hlahd.sh
RUN chmod +x /usr/local/bin/run_hlahd.sh

WORKDIR /data

ENTRYPOINT ["run_hlahd.sh"]
