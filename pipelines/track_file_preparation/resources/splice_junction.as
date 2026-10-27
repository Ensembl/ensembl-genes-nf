table spliceJunction
"Splice-junction intervals and read support. Coordinates use 0-based, half-open BED intervals."
(
    string chrom; "Reference sequence name"
    uint chromStart; "Junction start, 0-based"
    uint chromEnd; "Junction end, exclusive"
    string name; "Junction identifier"
    uint score; "Display score, 0-1000; capped from unique_reads"
    char[1] strand; "+, -, or ."
    uint unique_reads; "Number of uniquely mapping reads spanning the junction"
    uint multi_reads; "Number of multimapping reads spanning the junction"
    string splice_motif; "Splice-site dinucleotide pair, e.g. GT-AG, GC-AG, AT-AC; . means unavailable or non-canonical in STAR"
)
