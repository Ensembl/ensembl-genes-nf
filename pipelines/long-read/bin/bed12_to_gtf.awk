BEGIN { FS = "\t"; OFS = "\t" }
{
    if ($0 ~ /^#/ || NF == 0) next
    if (NF != 12) { print "invalid BED12 column count" > "/dev/stderr"; exit 1 }
    chrom = $1; start = $2; end = $3; tid = $4; strand = $6
    block_count = $10
    n_sizes = split($11, sizes, ",")
    n_offsets = split($12, offsets, ",")
    if (sizes[n_sizes] == "") n_sizes--
    if (offsets[n_offsets] == "") n_offsets--
    if (n_sizes != block_count || n_offsets != block_count || strand !~ /^[+-.]$/) {
        print "invalid BED12 block list" > "/dev/stderr"; exit 1
    }
    for (i = 1; i <= block_count; i++) {
        exon_start = start + offsets[i] + 1
        exon_end = exon_start + sizes[i] - 1
        if (sizes[i] <= 0 || exon_start < start + 1 || exon_end > end) {
            print "BED12 block outside transcript bounds" > "/dev/stderr"; exit 1
        }
        print chrom, "tama", "exon", exon_start, exon_end, ".", strand, ".", "gene_id \"" tid "\"; transcript_id \"" tid "\"; exon_number \"" i "\";"
    }
}
