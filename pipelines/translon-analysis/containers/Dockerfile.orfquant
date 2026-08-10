FROM quay.io/biocontainers/orfquant:1.1.0--r44h9ee0642_6

# ORFquant 1.1.0 predates the Bioconductor split of makeTxDbFromGFF into
# txdbmaker. Install the missing runtime package without changing ORFquant.
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_version("dbplyr", version="2.5.0", repos="https://cloud.r-project.org", upgrade="never"); BiocManager::install(c("BiocFileCache", "biomaRt", "txdbmaker"), ask=FALSE, update=FALSE); stopifnot(requireNamespace("txdbmaker", quietly=TRUE))' \
    && command -v ps >/dev/null
