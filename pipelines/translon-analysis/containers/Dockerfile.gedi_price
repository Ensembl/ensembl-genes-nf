FROM community.wave.seqera.io/library/gedi_price:2392624d5f803049

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends samtools \
    && rm -rf /var/lib/apt/lists/*
