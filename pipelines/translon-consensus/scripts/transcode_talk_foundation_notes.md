# Translon Consensus Talk Foundation Notes

## Talk spine

1. Callers produce different candidate universes.
   - The first plot uses the pooled pancreas context to compare like-for-like output scale and composition.
   - High call count is not treated as better; it defines each caller's search space.

2. Agreement is boundary-specific.
   - Full ORF agreement, shared start agreement, and shared stop agreement are different biological claims.
   - The notebook calculates all three in genomic coordinate space.

3. Identity is genomic; interpretation is transcript-relative.
   - Candidate identity is based on genomic coordinates and strand.
   - Context labels such as annotated CDS, non-CDS AUG, near-cognate, and mixed class are added afterward.

4. Disagreement has structure.
   - Boundary agreement and tool-specific calls are separated from full consensus.
   - Start-choice ambiguity is summarized through starts-per-stop tables.

5. Consensus must be tested against Ribo-seq evidence.
   - Consensus is not ground truth by itself.
   - Registered periodicity, coverage pass, uniformity, start rise, and stop drop-off are available in `evidence_df`.

6. P-site offset correction is required before interpreting periodicity.
   - The pancreas tracks show dominant CDS body frame 2 with fraction around 0.81-0.85.
   - The frozen offset is therefore 1 for all pancreas samples.
   - Old periodicity plots are invalid if they were generated before this offset was applied.

## Data cuts

- Main caller universe and agreement plots use `Ribo_Pancreas_pooled` where available.
- Ribo-seq evidence plots use pancreas scores and focus on `Ribo_Pancreas_pooled` for clean slide interpretation.
- The base DB filter is:
  - `qc_status = 'pass'`
  - `sample_id NOT GLOB '*_fastq'`
  - `source_feature_class IN ('cds', 'non_cds')`

## Notebook dataframes

- `calls_raw`: QC-pass DB rows.
- `calls_pooled`: deduplicated pooled `(tool, feature_key)` rows.
- `candidate_df`: one row per pooled genomic candidate.
- `agreement_df`: candidate-level full/start/stop support and agreement class.
- `context_df`: simple candidate context class.
- `scores_unique`: registered bigWig score table if available.
- `evidence_df`: merged candidate, agreement, context, and Ribo-seq evidence table.

## Output figures

- `01_caller_universe.png`: caller output scale and composition.
- `02_agreement_geometry.png`: full/start/stop agreement support tiers.
- `03_context_by_agreement.png`: context composition by agreement class.
- `04_starts_per_stop.png`: start-choice ambiguity at shared stops.
- `05_psite_offset_calibration.png`: P-site offset consistency across pancreas tracks.
- `06_registered_periodicity_by_class.png`: CDS vs non-CDS registered periodicity.
- `07_evidence_by_agreement.png`: registered periodicity by agreement class.
- `08_coverage_pass_by_tool.png`: fraction of calls with measurable Ribo-seq coverage.
- `09_locus_shared_stop_multiple_starts.png`: simple locus example for shared-stop/start-choice ambiguity.
- `10_locus_tool_specific_high_evidence.png`: simple locus example for a single-tool high-evidence call.
- `11_locus_consensus_low_evidence.png`: simple locus example for consensus with weak evidence.

## Caveats

- BigWig periodicity scores are invalid unless `signal_score_summary.json` records `psite_offset` for each scored sample.
- Do not reuse old `periodicity_distribution.png`, `12_evidence_vs_consensus_tier.png`, `13_consensus_x_evidence.png`, or `17_coverage_and_evidence_by_tool.png` if they came from unregistered scores.
- Start and stop flank metrics are currently genomic-flank based because the translon DB stores ORF blocks, not full transcript models for splice-aware UTR flanks.
- The simple context classes are intentionally coarse for tomorrow's talk; they are not a full transcript-relative biotype ontology.

## Commands to regenerate on `/hps`

```bash
cd /hps/software/users/ensembl/genebuild/jackt/ensembl-genes-nf/pipelines/translon-consensus/scripts
python3 -m jupyter nbconvert --to notebook --execute --inplace transcode_talk_foundation.ipynb
```

If the setup cell prints `Scorer supports P-site offsets: False`, sync the updated scorer first.
