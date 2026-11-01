# Run-level track schemas

This is the contract to review before processing the GDM backlog. The first
phase produces independent run-level tracks; it does not combine BAMs or
aggregate tracks. The downstream aggregation pipeline may consume these run
tracks later.

## Track artifacts

| Track type | Artifact | Coordinate/schema contract | Producer |
|---|---|---|---|
| `coverage_forward` | BigWig | 25 bp bins by default; deepTools `bamCoverage`, CPM-normalized by default, using `--filterRNAstrand forward` | STAR BAM + current BAM index |
| `coverage_reverse` | BigWig | Same as forward, using `--filterRNAstrand reverse` | STAR BAM + current BAM index |
| `gene_model` | BigBed | UCSC BED12, sorted by chromosome/start/end, created with `bedToBigBed -type=bed12 -tab` | manifest GTF or run-level Scallop GTF |
| `splice_junction` | BigBed | Nine-column `spliceJunction` AutoSql schema in `resources/splice_junction.as` | STAR `SJ.out.tab` |

The BED12 columns are `chrom`, `chromStart`, `chromEnd`, `name`, `score`,
`strand`, `thickStart`, `thickEnd`, `itemRgb`, `blockCount`, `blockSizes`, and
`blockStarts`. The gene-model path is identical whether its GTF came from the
manifest or from Scallop.

The splice-junction columns are `chrom`, `chromStart`, `chromEnd`, `name`,
`score`, `strand`, `unique_reads`, `multi_reads`, and `splice_motif`. STAR
coordinates are converted from 1-based inclusive starts to 0-based half-open
BED coordinates. The full field definitions are in `resources/splice_junction.as`.

## Scallop decision

For GDM run-level gene models, use the pinned Scallop container and leave the
library type at the provisional default (`scallop_library_type=auto`) until
the library protocol is confirmed. If the BAMs do not carry usable `XS` tags,
Scallop requires an explicit `unstranded`, `first`, or `second` setting. Scallop
emits one GTF per run; that GTF then passes through the existing validation and
GTF-to-BigBed chain.

The Scallop output is not a new public track format. It is the source for the
same `gene_model` BED12 BigBed and the same track-result/provenance schema.

## Required pre-processing decisions

Before launching the backlog, confirm these values for the GDM data:

1. The RNA-seq library type. `auto` is the safe operational default when the
   protocol metadata is not authoritative, but it should be replaced with
   `first`, `second`, or `unstranded` when known.
2. Whether the existing coverage convention is appropriate. It currently
   emits forward/reverse BigWigs using the dUTP/NSR/NNSR interpretation and
   CPM normalization with 25 bp bins.
3. The assembly-specific `chrom.sizes` files. Each per-GCA manifest must be
   run with the matching file; no cross-assembly size file should be reused.

The GDM manifest generator therefore defaults to run-level coverage and STAR
junction inputs for backwards compatibility. To prepare the complete set
including Scallop gene models, explicitly request:

```bash
--gene-model-source scallop \
--track-types coverage,gene_model,splice_junction
```

