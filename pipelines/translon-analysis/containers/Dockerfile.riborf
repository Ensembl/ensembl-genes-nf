FROM ubuntu:24.04
ARG RIBORF_COMMIT=be74bf9bf0b42b653d969a51f046cd38832136d9
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates git perl r-base samtools bowtie2 \
    && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/zhejilab/RibORF.git /opt/riborf \
    && cd /opt/riborf && git checkout ${RIBORF_COMMIT}
ENV RIBORF_HOME=/opt/riborf/RibORF.2.0
