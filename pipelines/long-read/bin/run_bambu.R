#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: run_bambu.R <input-bam> <reference-fasta> <output-gtf> <output-rds> <output-assignments.tsv> [track-reads] [quantify]")
suppressPackageStartupMessages(library(bambu))
bam <- normalizePath(args[[1]], mustWork = TRUE)
reference <- normalizePath(args[[2]], mustWork = TRUE)
output_gtf <- args[[3]]
output_rds <- args[[4]]
output_assignments <- args[[5]]
output_manifest <- if (length(args) >= 8) args[[8]] else "input_manifest.json"
track_reads <- if (length(args) >= 6) as.logical(args[[6]]) else FALSE
quantify <- if (length(args) >= 7) as.logical(args[[7]]) else TRUE

if (anyNA(c(track_reads, quantify))) stop("track-reads and quantify must be TRUE or FALSE")

sha256 <- function(path) {
    fields <- strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]]
    fields[[1]]
}
write_manifest <- function() {
    jsonlite::write_json(
        list(
            backend = "bambu",
            bam = list(path = bam, sha256 = sha256(bam)),
            reference = list(path = reference, sha256 = sha256(reference)),
            track_reads = track_reads,
            quantify = quantify
        ),
        output_manifest,
        auto_unbox = TRUE,
        pretty = TRUE
    )
}
genome <- Biostrings::readDNAStringSet(reference)
genome_names <- sub("[[:space:]].*$", "", names(genome))
if (anyDuplicated(genome_names)) stop("reference FASTA has duplicate sequence IDs after header normalization")
names(genome) <- genome_names
empty_result <- FALSE
result <- tryCatch(
    bambu::bambu(
        reads = bam,
        genome = genome,
        opt.discovery = list(fitReadClassModel = FALSE),
        discovery = TRUE,
        quant = quantify,
        trackReads = track_reads,
        ncore = 1,
        verbose = TRUE
    ),
    error = function(condition) {
        if (grepl("No valid annotations|No annotations were provided", conditionMessage(condition))) {
            empty_result <<- TRUE
            return(NULL)
        }
        stop(condition)
    }
)
if (empty_result) {
    writeLines("# EMPTY: Bambu produced no transcript models", output_gtf)
    saveRDS(NULL, output_rds)
    writeLines("read\\ttranscript", output_assignments)
    write_manifest()
    quit(save = "no", status = 0)
}

saveRDS(result, output_rds)
if (quantify) {
    output_dir <- file.path(getwd(), "bambu_output")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    bambu::writeBambuOutput(result, path = output_dir, prefix = "bambu")
    gtf_files <- list.files(output_dir, pattern = "\\.gtf$", full.names = TRUE)
    if (length(gtf_files) != 1) stop("Bambu did not produce exactly one GTF output")
    file.copy(gtf_files[[1]], output_gtf, overwrite = TRUE)
} else {
    bambu::writeToGTF(result, file = output_gtf)
}

assignments <- NULL
if (quantify && track_reads && "readToTranscriptMaps" %in% names(S4Vectors::metadata(result)))
    assignments <- S4Vectors::metadata(result)$readToTranscriptMaps
if (is.null(assignments)) {
    writeLines("read\\ttranscript", output_assignments)
} else {
    if (is.list(assignments)) assignments <- do.call(rbind, lapply(assignments, as.data.frame))
    utils::write.table(assignments, output_assignments, sep = "\\t", quote = FALSE, row.names = FALSE)
}
write_manifest()
