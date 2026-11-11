#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(CAGEfightR)
    library(GenomicRanges)
    library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
    index <- match(name, args)
    if (is.na(index) || index == length(args)) stop(sprintf("Missing argument %s", name))
    args[[index + 1L]]
}
get_args_after <- function(name) {
    index <- match(name, args)
    if (is.na(index) || index == length(args)) stop(sprintf("Missing argument %s", name))
    following <- args[(index + 1L):length(args)]
    next_option <- which(grepl("^--", following))[1L]
    if (is.na(next_option)) following else following[seq_len(next_option - 1L)]
}

assembly <- get_arg("--assembly")
species <- get_arg("--species")
metadata_path <- get_arg("--metadata")
gtf_path <- get_arg("--gtf")
output_dir <- get_arg("--output")
bed_paths <- get_args_after("--beds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(metadata_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(metadata) != length(bed_paths)) stop("Metadata rows and CTSS files differ")

read_ctss <- function(path) {
    values <- read.delim(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    if (ncol(values) < 6) stop(sprintf("CTSS BED must have six columns: %s", path))
    values <- values[, 1:6]
    colnames(values) <- c("chrom", "start", "end", "sample", "count", "strand")
    values$count <- as.numeric(values$count)
    values$key <- paste(values$chrom, values$start, values$strand, sep = ":")
    values
}

ctss_tables <- lapply(bed_paths, read_ctss)
all_keys <- sort(unique(unlist(lapply(ctss_tables, function(table) table$key))))
counts <- matrix(0, nrow = length(all_keys), ncol = length(ctss_tables), dimnames = list(all_keys, metadata$sample_id))
for (index in seq_along(ctss_tables)) {
    table <- ctss_tables[[index]]
    counts[match(table$key, all_keys), index] <- table$count
}
key_parts <- do.call(rbind, strsplit(all_keys, ":", fixed = TRUE))
row_ranges <- GRanges(seqnames = key_parts[, 1], ranges = IRanges(start = as.integer(key_parts[, 2]) + 1L, width = 1L), strand = key_parts[, 3])
ctss <- SummarizedExperiment(assays = list(counts = counts), rowRanges = row_ranges, colData = S4Vectors::DataFrame(metadata))
ctss <- calcTPM(ctss)
ctss <- calcPooled(ctss)
tss <- quickTSSs(ctss)
enhancers <- quickEnhancers(ctss)
tss <- quantifyClusters(ctss, rowRanges(tss), inputAssay = "counts", sparse = FALSE)
enhancers <- quantifyClusters(ctss, rowRanges(enhancers), inputAssay = "counts", sparse = FALSE)

write_bed <- function(ranges, path) {
    if (length(ranges) == 0L) {
        file.create(path)
        return(invisible(NULL))
    }
    output <- data.frame(as.character(seqnames(ranges)), start(ranges) - 1L, end(ranges), ifelse(is.na(ranges$cluster), "cluster", as.character(ranges$cluster)), ifelse(is.na(ranges$score), 0, as.numeric(ranges$score)), as.character(strand(ranges)))
    write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_bed(rowRanges(tss), file.path(output_dir, "tss_clusters.bed"))
write_bed(rowRanges(enhancers), file.path(output_dir, "enhancers.bed"))

annotation <- rtracklayer::import(gtf_path)
annotation <- annotation[annotation$type %in% c("gene", "transcript")]
hits <- findOverlaps(rowRanges(tss), annotation, ignore.strand = TRUE)
annotation_table <- data.frame(cluster_index = integer(), feature_type = character(), gene_id = character(), transcript_id = character(), stringsAsFactors = FALSE)
if (length(hits) > 0L) {
    query <- queryHits(hits)
    subject <- subjectHits(hits)
    get_attribute <- function(name, index) {
        values <- mcols(annotation)[[name]]
        if (is.null(values) || is.na(values[index])) "" else as.character(values[index])
    }
    annotation_table <- data.frame(
        cluster_index = query,
        feature_type = as.character(annotation$type[subject]),
        gene_id = vapply(subject, function(index) get_attribute("gene_id", index), character(1)),
        transcript_id = vapply(subject, function(index) get_attribute("transcript_id", index), character(1)),
        stringsAsFactors = FALSE
    )
}
write.table(annotation_table, file.path(output_dir, "tss_annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(assay(ctss, "counts")), file.path(output_dir, "ctss_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(assay(tss, "counts")), file.path(output_dir, "tss_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(assay(tss, "counts")), file.path(output_dir, "promoter_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(metadata, file.path(output_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
qc <- data.frame(sample_id = metadata$sample_id, library_size = colSums(counts), stringsAsFactors = FALSE)
write.table(qc, file.path(output_dir, "atlas.qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
