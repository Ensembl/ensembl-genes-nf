# CHM13 annotation release

This pipeline compares and reconciles projected and manual GFF3 annotations for
the CHM13 release. It remains a separate entrypoint from the whole-genome
projection workflow: projection produces a reproducible baseline, while manual
priority is an auditable downstream decision.

Before comparison, the workflow sanitizes the GFFs for the Ensembl core loader:
HAVANA `gene_segment` records are rewritten as `transcript`, and projected
`gene_type`/`transcript_type` values are copied to loader-compatible `biotype`
attributes. IDs, Parent links, sources, and all other attributes are retained.
The `gff_sanitization.*.json` reports record-level changes.

## Input CSV

```csv
sample_id,projected_gff,manual_gff,decision_tsv,assembly_fasta,assembly_fai
chm13,results/chm13/projected.gff3,havana/non_filler.gff3,-,chm13.fa,chm13.fa.fai
manual_only,-,havana/non_filler.gff3,-,chm13.fa,chm13.fa.fai
projected_only,results/chm13/projected.gff3,-,-,chm13.fa,chm13.fa.fai
```

Use `-` for an absent annotation source. The output contains:

- `annotation_review.tsv`: one row per manual transcript with structural class,
  overlapping projected models, a best-match explanation (`exons=...`,
  `five_prime_UTR=...`, `three_prime_UTR=...`), and the default action;
- `integrated.gff3`: projected models plus every manual model by default. Manual
  records retain their original IDs and parentage and have
  `source=HAVANA_manual`, `manual_source=HAVANA`,
  `manual_annotation_priority=gold`, `manual_review_required`,
  `manual_review_reason`, and `manual_recommended_action` attributes;
- `havana_decisions.tsv`: only the cases requiring HAVANA input, with blank
  `havana_decision` and `havana_comment` columns ready for handback;
- `review_cases.csv`: compact review coordinates for spreadsheet/JBrowse use;
- `review_cases_representative.csv`: a small balanced subset of novel, changed,
  UTR-different, and exact-match cases used by the generated session;
- `jbrowse_config.json` and `jbrowse_session.json`: assembly, tracks, and review
  views for a JBrowse 2 deployment.

For example, to explicitly promote one reviewed manual isoform:

```csv
manual_transcript_id,action
ENST_MANUAL_001,promote_canonical
```

The default policy adds every manual model unchanged alongside the projection.
It does not silently merge or discard manual records. Exact transcript matches
are marked as no-review-needed, while changed splice structures, ambiguous
parent-gene assignments, and readthrough transcripts are listed in
`havana_decisions.tsv`. To make a manual transcript replace projected
models, supply a decision table to the script with `keep_manual`; use
`promote_canonical` when the retained manual transcript should also be marked
`manual_canonical=true` and `manual_promotion_status=promoted`. Use
`keep_both`, `keep_projected`, or `skip` for explicit review decisions. The
decision-table hook is intentionally kept outside the automatic comparison so
gold manual annotation is never silently discarded by an inferred merge.

## Run

```bash
cd pipelines/chm13_annotation_release
nextflow run main.nf -profile slurm \\
  --sources_csv /path/sources.csv \\
  --assembly_name homo_sapiens_gca009914755v4 \\
  --outdir /path/results
```

The JBrowse JSON uses filenames relative to its output directory. Keep or copy
the assembly FASTA/FAI and the annotation GFF3 files beside the generated JSON
files, or edit the URIs in `jbrowse_config.json`. Serve the directory with a web
server or copy it into a JBrowse Web data directory; local `file://` URLs are
normally blocked by browsers.
