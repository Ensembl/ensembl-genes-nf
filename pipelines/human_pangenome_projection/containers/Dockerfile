# Container for the `hpp` pangenome projection tool.
#
# Bundles the hpp Python package plus the aligners it shells out to internally
# (gap-filling and per-gene rescue use minimap2; FASTA indexing uses samtools).
# The MAIN whole-genome alignment is NOT run here — that is a separate, swappable
# Nextflow module (MINIMAP2_WGA). This image is only for the projection stages.
#
# Build + push (then set params.hpp_container to the pushed URL):
#   docker build -t ghcr.io/ensembl/hpp:0.1.0 -f containers/Dockerfile .
#   docker push ghcr.io/ensembl/hpp:0.1.0

FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="hpp" \
      org.opencontainers.image.description="Human pangenome annotation projection tool" \
      org.opencontainers.image.version="0.1.0"

USER root

# System aligners / indexers used internally by the projection stages.
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        python=3.11 \
        minimap2=2.28 \
        samtools=1.20 \
        pysam \
        edlib \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# Install the hpp package (build context is the pipeline directory).
COPY pyproject.toml /opt/hpp/pyproject.toml
COPY hpp /opt/hpp/hpp
RUN pip install --no-cache-dir /opt/hpp

# Sanity check at build time.
RUN hpp --version

WORKDIR /work
