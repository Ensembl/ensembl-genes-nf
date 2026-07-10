FROM rocker/r-ver:4.3.3

ARG CHOROS_COMMIT=23d3e424f782d7b729d0c07989096cd00d6fa486

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git samtools python3 python3-pysam \
        libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev libxml2-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("BiocManager", "data.table", "doParallel", "foreach", "ggplot2", "MASS", "patchwork", "RColorBrewer", "reshape"), repos="https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install(c("Biostrings", "Rsamtools"), ask=FALSE, update=FALSE)'

# ChOROS declares vignette tooling as a hard runtime dependency upstream. The
# pipeline does not build vignettes, so remove only those two packaging-only
# dependencies from the pinned source before installation.
RUN git clone https://github.com/lareaulab/choros.git /tmp/choros \
    && cd /tmp/choros \
    && git checkout "${CHOROS_COMMIT}" \
    && sed -i 's/, knitr, rmarkdown, Biostrings/, Biostrings/' DESCRIPTION \
    && sed -i '/^VignetteBuilder:/d' DESCRIPTION \
    && R CMD INSTALL . \
    && cd / \
    && rm -rf /tmp/choros \
    && Rscript -e 'stopifnot(requireNamespace("choros", quietly=TRUE), requireNamespace("MASS", quietly=TRUE))'

LABEL org.opencontainers.image.source="https://github.com/lareaulab/choros"
LABEL org.opencontainers.image.description="Pinned ChOROS pilot runtime for the Ensembl Ribo-seq pipeline"
