#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Merge run-level GENCODE pancreas pilot bedGraphs into the six canonical pilot samples.

The default grouping is:
  SRR11005875..SRR11005879 -> SRR11005875_to_79
  SRR11005880..SRR11005884 -> SRR11005880_to_84
  SRR11005885..SRR11005889 -> SRR11005885_to_89
  SRR11005890..SRR11005894 -> SRR11005890_to_94
  SRR11005895..SRR11005899 -> SRR11005895_to_99
  SRR11005900..SRR11005904 -> SRR11005900_to_04

Usage:
  merge_gencode_pilot_bedgraphs.sh --bedgraphs DIR --outdir DIR [options]

Required:
  --bedgraphs DIR       Directory containing run bedGraphs from riboseq output.
  --outdir DIR          Output directory for grouped bedGraphs and optional BigWigs.

Options:
  --chrom-sizes FILE    Convert grouped bedGraphs to BigWig using this chrom sizes file.
  --manifest FILE       Two-column TSV: run<TAB>pilot_sample. Defaults to built-in pancreas grouping.
  --types LIST          Comma-separated BAM types to merge.
                        Default: unique_no_junction,unique_with_junction,multi_no_junction,multi_with_junction
  --strands LIST        Comma-separated strands to merge. Default: forward,reverse
  --aggregation MODE    sum or mean. Default: sum. Use sum for pooled-signal pilot tracks.
  --force               Overwrite existing outputs.
  --dry-run             Print planned operations without writing outputs.
  -h, --help            Show this help.

Outputs:
  OUTDIR/bedgraphs/<pilot_sample>.<bam_type>.<strand>.sorted.bedgraph
  OUTDIR/bigwigs/<pilot_sample>.<bam_type>.<strand>.bw  (if --chrom-sizes is provided)
  OUTDIR/pilot_group_manifest.tsv
EOF
}

bedgraphs_dir=""
outdir=""
chrom_sizes=""
manifest=""
types_csv="unique_no_junction,unique_with_junction,multi_no_junction,multi_with_junction"
strands_csv="forward,reverse"
aggregation="sum"
force=0
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bedgraphs)
      bedgraphs_dir="${2:?missing value for --bedgraphs}"
      shift 2
      ;;
    --outdir)
      outdir="${2:?missing value for --outdir}"
      shift 2
      ;;
    --chrom-sizes)
      chrom_sizes="${2:?missing value for --chrom-sizes}"
      shift 2
      ;;
    --manifest)
      manifest="${2:?missing value for --manifest}"
      shift 2
      ;;
    --types)
      types_csv="${2:?missing value for --types}"
      shift 2
      ;;
    --strands)
      strands_csv="${2:?missing value for --strands}"
      shift 2
      ;;
    --aggregation)
      aggregation="${2:?missing value for --aggregation}"
      shift 2
      ;;
    --force)
      force=1
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$bedgraphs_dir" || -z "$outdir" ]]; then
  echo "ERROR: --bedgraphs and --outdir are required" >&2
  usage >&2
  exit 2
fi

if [[ "$aggregation" != "sum" && "$aggregation" != "mean" ]]; then
  echo "ERROR: --aggregation must be 'sum' or 'mean'" >&2
  exit 2
fi

if [[ ! -d "$bedgraphs_dir" ]]; then
  echo "ERROR: bedGraph directory does not exist: $bedgraphs_dir" >&2
  exit 1
fi

if [[ -n "$chrom_sizes" && ! -s "$chrom_sizes" ]]; then
  echo "ERROR: chrom sizes file does not exist or is empty: $chrom_sizes" >&2
  exit 1
fi

if [[ "$dry_run" -eq 0 ]] && ! command -v bedtools >/dev/null 2>&1; then
  echo "ERROR: bedtools is required for bedtools unionbedg" >&2
  exit 1
fi

if [[ "$dry_run" -eq 0 && -n "$chrom_sizes" ]] && ! command -v bedGraphToBigWig >/dev/null 2>&1; then
  echo "ERROR: bedGraphToBigWig is required when --chrom-sizes is provided" >&2
  exit 1
fi

IFS=',' read -r -a bam_types <<< "$types_csv"
IFS=',' read -r -a strands <<< "$strands_csv"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

manifest_work="$tmpdir/pilot_group_manifest.tsv"
if [[ -n "$manifest" ]]; then
  if [[ ! -s "$manifest" ]]; then
    echo "ERROR: manifest does not exist or is empty: $manifest" >&2
    exit 1
  fi
  awk 'BEGIN{FS=OFS="\t"} NR==1 && $1=="run" {next} NF>=2 {print $1,$2}' "$manifest" > "$manifest_work"
else
  cat > "$manifest_work" <<'EOF'
SRR11005875	SRR11005875_to_79
SRR11005876	SRR11005875_to_79
SRR11005877	SRR11005875_to_79
SRR11005878	SRR11005875_to_79
SRR11005879	SRR11005875_to_79
SRR11005880	SRR11005880_to_84
SRR11005881	SRR11005880_to_84
SRR11005882	SRR11005880_to_84
SRR11005883	SRR11005880_to_84
SRR11005884	SRR11005880_to_84
SRR11005885	SRR11005885_to_89
SRR11005886	SRR11005885_to_89
SRR11005887	SRR11005885_to_89
SRR11005888	SRR11005885_to_89
SRR11005889	SRR11005885_to_89
SRR11005890	SRR11005890_to_94
SRR11005891	SRR11005890_to_94
SRR11005892	SRR11005890_to_94
SRR11005893	SRR11005890_to_94
SRR11005894	SRR11005890_to_94
SRR11005895	SRR11005895_to_99
SRR11005896	SRR11005895_to_99
SRR11005897	SRR11005895_to_99
SRR11005898	SRR11005895_to_99
SRR11005899	SRR11005895_to_99
SRR11005900	SRR11005900_to_04
SRR11005901	SRR11005900_to_04
SRR11005902	SRR11005900_to_04
SRR11005903	SRR11005900_to_04
SRR11005904	SRR11005900_to_04
EOF
fi

if [[ "$dry_run" -eq 0 ]]; then
  mkdir -p "$outdir/bedgraphs"
  if [[ -n "$chrom_sizes" ]]; then
    mkdir -p "$outdir/bigwigs"
  fi
  {
    printf 'run\tpilot_sample\n'
    cat "$manifest_work"
  } > "$outdir/pilot_group_manifest.tsv"
fi

pilot_samples=()
while IFS= read -r pilot_sample; do
  pilot_samples+=("$pilot_sample")
done < <(awk -F '\t' '{print $2}' "$manifest_work" | sort -u)

missing_count=0
merged_count=0
bigwig_count=0

for pilot_sample in "${pilot_samples[@]}"; do
  runs=()
  while IFS= read -r run; do
    runs+=("$run")
  done < <(awk -F '\t' -v sample="$pilot_sample" '$2 == sample {print $1}' "$manifest_work")
  if [[ "${#runs[@]}" -eq 0 ]]; then
    continue
  fi

  for bam_type in "${bam_types[@]}"; do
    for strand in "${strands[@]}"; do
      inputs=()
      for run in "${runs[@]}"; do
        bg="$bedgraphs_dir/${run}.${bam_type}.${strand}.sorted.bedgraph"
        if [[ -s "$bg" ]]; then
          inputs+=("$bg")
        else
          echo "WARN: missing or empty input: $bg" >&2
          missing_count=$((missing_count + 1))
        fi
      done

      if [[ "${#inputs[@]}" -eq 0 ]]; then
        echo "WARN: no inputs for ${pilot_sample}.${bam_type}.${strand}; skipping" >&2
        continue
      fi

      out_bg="$outdir/bedgraphs/${pilot_sample}.${bam_type}.${strand}.sorted.bedgraph"
      out_bw="$outdir/bigwigs/${pilot_sample}.${bam_type}.${strand}.bw"

      if [[ "$force" -eq 0 && -e "$out_bg" ]]; then
        echo "ERROR: output exists, use --force to overwrite: $out_bg" >&2
        exit 1
      fi
      if [[ -n "$chrom_sizes" && "$force" -eq 0 && -e "$out_bw" ]]; then
        echo "ERROR: output exists, use --force to overwrite: $out_bw" >&2
        exit 1
      fi

      echo "Merging ${#inputs[@]} bedGraphs -> $out_bg"
      if [[ "$dry_run" -eq 0 ]]; then
        if [[ "$aggregation" == "sum" ]]; then
          bedtools unionbedg -filler 0 -i "${inputs[@]}" \
            | awk 'BEGIN{OFS="\t"} {sum=0; for (i=4; i<=NF; i++) sum += $i; print $1,$2,$3,sum}' \
            | sort -k1,1 -k2,2n > "$out_bg"
        else
          bedtools unionbedg -filler 0 -i "${inputs[@]}" \
            | awk 'BEGIN{OFS="\t"} {sum=0; n=NF-3; for (i=4; i<=NF; i++) sum += $i; print $1,$2,$3,(n ? sum/n : 0)}' \
            | sort -k1,1 -k2,2n > "$out_bg"
        fi
      fi
      merged_count=$((merged_count + 1))

      if [[ -n "$chrom_sizes" ]]; then
        echo "Converting -> $out_bw"
        if [[ "$dry_run" -eq 0 ]]; then
          bedGraphToBigWig "$out_bg" "$chrom_sizes" "$out_bw"
        fi
        bigwig_count=$((bigwig_count + 1))
      fi
    done
  done
done

echo "Done."
echo "Grouped bedGraphs planned/written: $merged_count"
echo "Grouped BigWigs planned/written: $bigwig_count"
echo "Missing input bedGraphs observed: $missing_count"
