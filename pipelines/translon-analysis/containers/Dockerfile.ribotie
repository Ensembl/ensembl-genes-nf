FROM python:3.11-slim
ARG RIBOTIE_COMMIT=0d26f2ed5c5642a8ab6a5d63c40360a2210ebd72
RUN apt-get update && apt-get install -y --no-install-recommends git gcc g++ \
    && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/TRISTAN-ORF/RiboTIE.git /opt/ribotie \
    && cd /opt/ribotie && git checkout ${RIBOTIE_COMMIT} \
    && pip install --no-cache-dir --index-url https://download.pytorch.org/whl/cpu torch==2.4.1 \
    && pip install --no-cache-dir transcript-transformer
ENV PATH=/opt/ribotie:$PATH
