# Pilot diagnostics

`diagnose_pilot.py` joins the published RDP selection, getRPF extraction
reports, QC metrics, gate decisions, and sample-completeness table.  It is
read-only and writes one row per run, with a diagnosis such as
`withheld_at_structure_inference`, `emitted_but_qc_gate_failed`,
`emitted_but_qc_missing_or_software_failure`, or `pass_good_gate`.

Run it against a published pilot with absolute paths:

```bash
python3 \
  /hps/software/users/ensembl/genebuild/jackt/ensembl-genes-nf/pipelines/riboseq/bin/diagnose_pilot.py \
  --pilot-root /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/rdp-production/release_2026-11-v4/pilot/danio_rerio \
  --selection /homes/jackt/rdp-production/inputs/release_2026-11-v4/pilot/selection.tsv \
  --output /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/rdp-production/release_2026-11-v4/pilot/danio_rerio/pipeline_info/pilot_diagnostics.tsv
```

The adjacent `.summary.json` reports diagnosis counts.  The report does not
reinterpret QC thresholds or alter read trimming; it records the evidence
needed to decide whether a protocol-specific reprocessing experiment is
warranted.
