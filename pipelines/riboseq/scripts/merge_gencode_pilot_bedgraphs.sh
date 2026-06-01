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
                        Alias: --input-dir.
  --outdir DIR          Output directory for grouped bedGraphs and optional BigWigs.
                        Alias: --output-dir.

Options:
  --chrom-sizes FILE    Convert grouped bedGraphs to BigWig using this chrom sizes file.
  --manifest FILE       Two-column TSV: run<TAB>pilot_sample. Defaults to built-in pancreas grouping.
  --skip-runs LIST      Comma-separated runs to skip cleanly. Missing or zero-byte bedGraphs
                        for these runs are not counted as incomplete inputs.
  --types LIST          Comma-separated BAM types to merge.
                        Default: unique_no_junction,unique_with_junction,multi_no_junction,multi_with_junction
  --strands LIST        Comma-separated strands to merge. Default: forward,reverse
  --aggregation MODE    sum or mean. Default: sum. Use sum for pooled-signal pilot tracks.
  --require-complete-groups
                        Fail if any non-skipped group/type/strand input file is missing or empty.
  --force               Overwrite existing outputs.
  --validate-existing   Check existing grouped outputs and report whether they are non-empty or valid empty files.
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
skip_runs_csv=""
types_csv="unique_no_junction,unique_with_junction,multi_no_junction,multi_with_junction"
strands_csv="forward,reverse"
aggregation="sum"
force=0
dry_run=0
require_complete_groups=0
validate_existing=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bedgraphs|--input-dir)
      bedgraphs_dir="${2:?missing value for --bedgraphs}"
      shift 2
      ;;
    --outdir|--output-dir)
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
    --skip-runs)
      skip_runs_csv="${2:?missing value for --skip-runs}"
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
    --validate-existing)
      validate_existing=1
      shift
      ;;
    --require-complete-groups)
      require_complete_groups=1
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
  echo "ERROR: --bedgraphs/--input-dir and --outdir/--output-dir are required" >&2
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

if [[ "$dry_run" -eq 0 && -n "$chrom_sizes" ]] && ! command -v bedGraphToBigWig >/dev/null 2>&1; then
  echo "ERROR: bedGraphToBigWig is required when --chrom-sizes is provided" >&2
  exit 1
fi

IFS=',' read -r -a bam_types <<< "$types_csv"
IFS=',' read -r -a strands <<< "$strands_csv"
IFS=',' read -r -a skip_runs <<< "$skip_runs_csv"

is_skip_run() {
  local query="$1"
  local skip_run
  for skip_run in "${skip_runs[@]}"; do
    if [[ -n "$skip_run" && "$query" == "$skip_run" ]]; then
      return 0
    fi
  done
  return 1
}

join_by_comma() {
  local IFS=','
  echo "$*"
}

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
empty_input_count=0
incomplete_count=0
skipped_run_count=0
valid_existing_output_count=0
empty_existing_output_count=0
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
      skipped=()
      expected_count=0
      missing_or_empty_count=0
      for run in "${runs[@]}"; do
        if is_skip_run "$run"; then
          skipped+=("$run")
          continue
        fi
        expected_count=$((expected_count + 1))
        bg="$bedgraphs_dir/${run}.${bam_type}.${strand}.sorted.bedgraph"
        if [[ -s "$bg" ]]; then
          inputs+=("$bg")
        elif [[ -e "$bg" ]]; then
          echo "WARN: zero-byte input excluded from merge: $bg" >&2
          empty_input_count=$((empty_input_count + 1))
          missing_or_empty_count=$((missing_or_empty_count + 1))
        else
          echo "WARN: missing input: $bg" >&2
          missing_count=$((missing_count + 1))
          missing_or_empty_count=$((missing_or_empty_count + 1))
        fi
      done

      if [[ "${#skipped[@]}" -gt 0 ]]; then
        echo "SKIP: ${pilot_sample}.${bam_type}.${strand}: skipped runs $(join_by_comma "${skipped[@]}")"
        skipped_run_count=$((skipped_run_count + ${#skipped[@]}))
      fi

      if [[ "$missing_or_empty_count" -gt 0 ]]; then
        incomplete_count=$((incomplete_count + 1))
        if [[ "$require_complete_groups" -eq 1 ]]; then
          echo "ERROR: incomplete group for ${pilot_sample}.${bam_type}.${strand}: ${#inputs[@]}/${expected_count} non-skipped inputs are non-empty" >&2
          exit 1
        fi
      fi

      out_bg="$outdir/bedgraphs/${pilot_sample}.${bam_type}.${strand}.sorted.bedgraph"
      out_bw="$outdir/bigwigs/${pilot_sample}.${bam_type}.${strand}.bw"

      if [[ "$validate_existing" -eq 1 ]]; then
        if [[ -s "$out_bg" ]]; then
          echo "VALID: existing non-empty grouped bedGraph: $out_bg"
          valid_existing_output_count=$((valid_existing_output_count + 1))
        elif [[ -e "$out_bg" ]]; then
          echo "VALID_EMPTY: existing empty grouped bedGraph: $out_bg"
          empty_existing_output_count=$((empty_existing_output_count + 1))
        fi
        if [[ -n "$chrom_sizes" && -s "$out_bw" ]]; then
          echo "VALID: existing BigWig: $out_bw"
          valid_existing_output_count=$((valid_existing_output_count + 1))
        fi
      fi

      if [[ "$dry_run" -eq 0 && "$force" -eq 0 && -e "$out_bg" ]]; then
        echo "ERROR: output exists, use --force to overwrite: $out_bg" >&2
        exit 1
      fi
      if [[ "$dry_run" -eq 0 && -n "$chrom_sizes" && "$force" -eq 0 && -e "$out_bw" ]]; then
        echo "ERROR: output exists, use --force to overwrite: $out_bw" >&2
        exit 1
      fi

      if [[ "${#inputs[@]}" -eq 0 ]]; then
        echo "WARN: no non-empty inputs for ${pilot_sample}.${bam_type}.${strand}; writing empty grouped bedGraph: $out_bg" >&2
      fi

      echo "Merging ${#inputs[@]} non-empty bedGraphs -> $out_bg"
      if [[ "$dry_run" -eq 0 ]]; then
        if [[ "${#inputs[@]}" -eq 0 ]]; then
          : > "$out_bg"
        else
          if ! command -v bedtools >/dev/null 2>&1; then
            echo "ERROR: bedtools is required for bedtools unionbedg" >&2
            exit 1
          fi
          if [[ "$aggregation" == "sum" ]]; then
            bedtools unionbedg -filler 0 -i "${inputs[@]}" \
              | awk 'BEGIN{OFS="\t"} {sum=0; for (i=4; i<=NF; i++) sum += $i; print $1,$2,$3,sum}' \
              | sort -k1,1 -k2,2n > "$out_bg"
          else
            bedtools unionbedg -filler 0 -i "${inputs[@]}" \
              | awk -v n_inputs="${#inputs[@]}" 'BEGIN{OFS="\t"} {sum=0; for (i=4; i<=NF; i++) sum += $i; print $1,$2,$3,(n_inputs ? sum/n_inputs : 0)}' \
              | sort -k1,1 -k2,2n > "$out_bg"
          fi
        fi
      fi
      merged_count=$((merged_count + 1))

      if [[ -n "$chrom_sizes" ]]; then
        if [[ "$dry_run" -eq 0 && ! -s "$out_bg" ]]; then
          echo "INFO: grouped bedGraph is empty; skipping BigWig conversion: $out_bg" >&2
        else
          echo "Converting -> $out_bw"
          if [[ "$dry_run" -eq 0 ]]; then
            bedGraphToBigWig "$out_bg" "$chrom_sizes" "$out_bw"
          fi
          bigwig_count=$((bigwig_count + 1))
        fi
      fi
    done
  done
done

echo "Done."
echo "Grouped bedGraphs planned/written: $merged_count"
echo "Grouped BigWigs planned/written: $bigwig_count"
echo "Missing input bedGraphs observed: $missing_count"
echo "Zero-byte non-skipped input bedGraphs observed: $empty_input_count"
echo "Skipped run/type/strand inputs observed: $skipped_run_count"
echo "Incomplete group/type/strand merges observed: $incomplete_count"
echo "Valid existing outputs observed: $valid_existing_output_count"
echo "Valid empty existing grouped bedGraphs observed: $empty_existing_output_count"
