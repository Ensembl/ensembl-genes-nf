BEGIN {
    OFS = "\t"
    source_alignments = 0
    inflated_alignments = 0
    suffixed_alignments = 0
    failed = 0
}

/^@/ {
    print
    next
}

{
    qname = $1
    copies = 1

    if (qname ~ /_x[0-9]+$/) {
        match(qname, /_x[0-9]+$/)
        copies = substr(qname, RSTART + 2) + 0
        if (copies < 1) {
            print "Invalid multiplicity in read name: " qname > "/dev/stderr"
            failed = 1
            exit 2
        }
        suffixed_alignments++
    } else if (qname ~ /_x/) {
        print "Malformed multiplicity suffix in read name: " qname > "/dev/stderr"
        failed = 1
        exit 2
    }

    source_alignments++
    for (i = 1; i <= copies; i++) {
        print
        inflated_alignments++
    }
}

END {
    if (!failed) {
        print "metric" OFS "value" > manifest
        print "source_alignments" OFS source_alignments >> manifest
        print "inflated_alignments" OFS inflated_alignments >> manifest
        print "suffixed_alignments" OFS suffixed_alignments >> manifest
        close(manifest)
    }
}
