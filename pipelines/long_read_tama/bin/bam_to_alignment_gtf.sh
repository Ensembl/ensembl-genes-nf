#!/usr/bin/env bash
set -euo pipefail

bam=${1:?BAM input is required}
output=${2:?GTF output is required}

# Keep this converter runnable in the samtools image, which intentionally does
# not include Python. samtools preserves coordinate order; sort the generated
# GTF explicitly before handing it to tmerge.
samtools view -F 2308 "$bam" \
    | awk '
        BEGIN { OFS = "\t" }
        function emit_exon(contig, exon_start, exon_end, strand, read_id, exon_number, attrs) {
            if (exon_start < exon_end) {
                attrs = "transcript_id \"" read_id "." seen[read_id] "\"; exon_number \"" exon_number "\"; read_id \"" read_id "\";"
                print contig, "tmerge", "exon", exon_start, exon_end, ".", strand, ".", attrs
            }
        }
        {
            read_id = $1
            seen[read_id]++
            contig = $3
            start = $4 + 0
            cigar = $6
            strand = (int($2 / 16) % 2) ? "-" : "+"
            ref = start
            exon_start = ref
            exon_number = 0

            while (cigar != "" && match(cigar, /[0-9]+[MIDNSHP=X]/)) {
                token = substr(cigar, RSTART, RLENGTH)
                length_bp = token + 0
                operation = substr(token, length(token), 1)
                cigar = substr(cigar, RSTART + RLENGTH)

                if (operation == "N") {
                    if (exon_start < ref) {
                        exon_number++
                        emit_exon(contig, exon_start, ref, strand, read_id, exon_number)
                    }
                    ref += length_bp
                    exon_start = ref
                } else if (operation == "M" || operation == "D" || operation == "=" || operation == "X") {
                    ref += length_bp
                } else if (operation == "S" || operation == "H" || operation == "P" || operation == "I") {
                    # These operations do not consume reference coordinates.
                }
            }
            if (exon_start < ref) {
                exon_number++
                emit_exon(contig, exon_start, ref, strand, read_id, exon_number)
            }
        }
    ' \
    | LC_ALL=C sort -k1,1 -k4,4n -k5,5n -k7,7 -k9,9 > "$output"

test -s "$output"
