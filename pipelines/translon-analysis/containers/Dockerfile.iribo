FROM ubuntu:24.04
ARG IRIBO_COMMIT=b2bfd3d85b64bbe770a045ee027065b22e0f5aaf
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ca-certificates git r-base r-cran-argparse r-cran-scales r-cran-ggplot2 r-cran-future.apply samtools \
    && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/CarvunisLab/iRibo.git /opt/iribo \
    && cd /opt/iribo && git checkout ${IRIBO_COMMIT} \
    && sed -i '20i#include <unordered_map>' iRibo.cpp \
    && make
ENV PATH=/opt/iribo:$PATH
ENV IRIBO_HOME=/opt/iribo
