#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(choros))
suppressPackageStartupMessages(library(MASS))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9) {
  stop("usage: run_choros.R BAM FASTA LENGTHS OFFSETS PREFIX NUM_GENES MIN_COVERAGE MIN_NONZERO CORES")
}

bam_fname <- args[[1]]
fasta_fname <- args[[2]]
lengths_fname <- args[[3]]
offsets_fname <- args[[4]]
prefix <- args[[5]]
num_genes <- as.integer(args[[6]])
min_coverage <- as.numeric(args[[7]])
min_nonzero <- as.integer(args[[8]])
num_cores <- as.integer(args[[9]])

alignment <- load_bam(
  bam_fname, fasta_fname, lengths_fname, offsets_fname,
  f5_length=3, f3_length=3, num_cores=num_cores
)
if (nrow(alignment) == 0) stop("no ChOROS-compatible alignments remain")

digest <- count_d5_d3(alignment)
subsets <- choose_subsets(digest, min_prop=0.9)[, c("d5", "d3")]
coverage <- calculate_transcript_density(alignment, lengths_fname)
nonzero <- count_nonzero_codons(alignment)
nonzero_for_coverage <- nonzero[names(coverage)]
nonzero_for_coverage[is.na(nonzero_for_coverage)] <- 0
eligible <- names(coverage)[
  coverage > min_coverage & nonzero_for_coverage > min_nonzero
]
training_set <- head(eligible, num_genes)
if (length(training_set) < 2) {
  stop("fewer than two transcripts satisfy ChOROS training thresholds")
}

training <- init_data(
  fasta_fname, lengths_fname, d5_d3_subsets=subsets,
  f5_length=3, f3_length=3, which_transcripts=training_set,
  num_cores=num_cores
)
training$transcript <- relevel(training$transcript, ref=training_set[[1]])
training$count <- count_footprints(alignment, training, "count")
model <- formula(count ~ transcript + A + P + E + d5*f5 + d3*f3 + gc)
fit <- MASS::glm.nb(model, data=training, model=FALSE)
coefs <- parse_coefs(fit)
if (any(!is.finite(coefs$estimate))) stop("non-finite ChOROS regression coefficient")

alignment$corrected <- correct_bias(alignment, coefs)
write.table(
  alignment, gzfile(paste0(prefix, ".choros_counts.tsv.gz")),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  coefs, paste0(prefix, ".choros_coefficients.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  data.frame(
    metric=c("alignments", "training_transcripts", "raw_count", "corrected_count"),
    value=c(nrow(alignment), length(training_set), sum(alignment$count), sum(alignment$corrected))
  ),
  paste0(prefix, ".choros_metrics.tsv"), sep="\t", quote=FALSE, row.names=FALSE
)
