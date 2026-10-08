# DIAMOND QC

Integrated runs use the validated backend model channels, so no post-run
annotation manifest is needed.

`--diamond_scope` controls the model scope:

- `cohort` (default): run one QC branch per finalized cohort/backend model.
- `accession`: run one QC branch per finalized accession/backend model.
- `both`: run both scopes.

Accession inputs are the same validated BED12 products used by accession model
comparison. The workflow adds explicit `scope` metadata before invoking the
shared `RUN_DIAMOND_QC` subworkflow, preserving accession, backend, and scope
in task tags and report names.

## QC stages

1. `EXTRACT_COMBINED_TRANSCRIPTS` runs `bedtools getfasta -split -s` to extract
   spliced transcript sequences from each BED12 model.
2. `PREDICT_LONGEST_ATG_ORFS` runs `longest_atg_orf.py` and emits a peptide FASTA
   plus an ORF manifest. Models without a usable peptide remain represented as
   `NO_ATG` or `PARTIAL_ORF` rows.
3. `DIAMOND_BLASTP` runs Diamond 2.1.24 with `--evalue 1e-5`,
   `--max-target-seqs 1`, and tabular output containing identity, alignment
   length, query coverage, e-value, and bit score.
4. `REPORT_COMBINED_DIAMOND_MODELS` left-joins the best hit to every ORF row and
   emits `CLASSIFIED`, `NO_PROTEIN_HIT`, or `NO_PEPTIDE` status. Classified hits
   also receive the selected legacy coverage/identity class.

This is protein-homology support for predicted ORFs, not transcript-level
functional annotation. A missing hit can reflect no complete ORF, no database
match, or limitations of the longest-ATG heuristic.

The existing `--diamond_annotation_manifest` mode remains available when
annotations already exist outside the integrated workflow.
