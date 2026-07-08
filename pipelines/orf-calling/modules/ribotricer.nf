/* ribotricer module: PREP -> RUN -> PARSE */

process PREP_RIBOTRICER {
  tag "${meta.id}"
  label 'process_low'
  conda "bioconda::ribotricer"
  container "${ params.container_ribotricer ?: 'quay.io/biocontainers/ribotricer:latest' }"
  publishDir "${params.outdir}/orf_calls/ribotricer/${meta.id}", mode: 'copy', pattern: "prep/**"

  input:
  tuple val(meta), path(bam), path(bai)
  path gtf
  path fasta

  output:
  tuple val(meta), path('prep'), emit: prepared

  script:
  """
  mkdir -p prep
  prep_inputs.py --outdir prep --sample-id ${meta.id} \
    --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} \
    --bam-type transcriptome
  """
  stub:
  """
  mkdir -p prep && echo '{}' > prep/inputs.json
  """
}

process RUN_RIBOTRICER {
  tag "${meta.id}"
  label 'process_medium'
  conda "bioconda::ribotricer"
  container "${ params.container_ribotricer ?: 'quay.io/biocontainers/ribotricer:latest' }"
  publishDir "${params.outdir}/orf_calls/ribotricer/${meta.id}", mode: 'copy', pattern: "raw/**"

  input:
  tuple val(meta), path(prepdir)

  output:
  tuple val(meta), path('raw'), emit: raw

  script:
  """
  mkdir -p raw
  export THREADS=${params.threads_ribotricer ?: 4}
  if [ -n "${task.stub}" ]; then
    echo -e "transcript_id\tstart\tend\tframe\tscore\tperiodicity" > raw/ribotricer.tsv
    echo -e "TX1\t120\t210\t0\t5.2\t0.8" >> raw/ribotricer.tsv
  else
    python3 - <<'PY'
import json, subprocess, os
from pathlib import Path
inp=json.loads(Path('prep/inputs.json').read_text())
bam,gtf,fasta = inp['bam'], inp['gtf'], inp['fasta']
threads=int(os.environ.get('THREADS','4'))
Path('raw').mkdir(exist_ok=True)
subprocess.run(['ribotricer','prepare-orfs','-a',gtf,'-g',fasta,'-o','raw/orfs'], check=False)
idx = 'raw/orfs_candidate_ORFs.tsv'
cmd=['ribotricer','detect-orfs','-b',bam,'-i',idx,'-o','raw','-t',str(threads)]
pl=inp.get('pass_lengths')
if pl and Path(pl).exists():
    L=[l.strip() for l in Path(pl).read_text().splitlines()[1:] if l.strip()]
    if L:
        cmd += ['--read_lengths', ','.join(L)]
subprocess.run(cmd, check=False)
PY
  fi
  """
  stub:
  """
  mkdir -p raw && touch raw/ribotricer.tsv
  """
}

process PARSE_RIBOTRICER {
  tag "${meta.id}"
  label 'process_low'
  container "ubuntu:22.04"
  publishDir "${params.outdir}/orf_calls/ribotricer/${meta.id}", mode: 'copy', pattern: "*"

  input:
  tuple val(meta), path(rawdir)

  output:
  tuple val(meta), path("${meta.id}.orf_calls.tsv"), emit: standardized
  tuple val(meta), path("${meta.id}.orf_calls.bed12"), emit: bed12
  path "versions.yml", emit: versions

  script:
  """
  python3 - << 'PY'
import csv, json
from pathlib import Path
from parse_common import write_header, append_rows, default_row
meta = json.loads('''${meta.toString().replace("'","\'")}''')
raw = Path('${rawdir}')
out = Path(f"{meta['id']}.orf_calls.tsv")
write_header(out)
fp = raw/'ribotricer.tsv'
rows = []
if fp.exists():
  r = csv.DictReader(open(fp), delimiter='\t')
  for row in r:
    rows.append([meta['id'],'ribotricer',None,row['start'],row['end'],None,row['frame'],row['transcript_id'],None,row['score'],None,None,json.dumps({'periodicity': row.get('periodicity')})])
else:
  rows.append(default_row(meta['id'],'ribotricer'))
append_rows(out, rows)
PY
  # Convert TSV to BED12 (minimal, single-block)
  python3 - << 'PY'
import csv, math
from pathlib import Path

def clamp(v, lo=0, hi=1000):
    return max(lo, min(hi, v))

def score_from_fields(score_str, pval_str):
    try:
        if pval_str and pval_str.lower() != 'none':
            p = float(pval_str)
            if math.isfinite(p) and 0 <= p <= 1:
                return clamp(int(round(1000 * (1.0 - p))))
    except Exception:
        pass
    try:
        if score_str and score_str.lower() != 'none':
            s = float(score_str)
            if math.isfinite(s):
                y = 1.0 - math.exp(-abs(s))
                return clamp(int(round(1000 * y)))
    except Exception:
        pass
    return 0

tsv = Path("${meta.id}.orf_calls.tsv")
bed = Path("${meta.id}.orf_calls.bed12")
with open(tsv, newline='') as fh, open(bed, 'w') as out:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        try:
            start = int(row['start']); end = int(row['end'])
        except Exception:
            continue
        if end <= start:
            continue
        chrom = (row.get('chrom') or '').strip() or (row.get('transcript_id') or 'unknown')
        name = f"{row.get('orf_id') or 'orf'}|{row.get('transcript_id') or 'tx'}|{row.get('tool') or 'tool'}"
        score = score_from_fields(row.get('score'), row.get('pval'))
        strand = (row.get('strand') or '.').strip() or '.'
        bed12 = '\t'.join([
            chrom, str(start), str(end), name, str(score), strand,
            str(start), str(end), '0,0,0', '1', f"{end-start},", '0,'
        ])
        out.write(bed12 + '\n')
PY
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      parser: python-3.11
  END_VERSIONS
  """
  stub:
  """
  echo -e "sample_id\ttool\tchrom\tstart\tend\tstrand\tframe\ttranscript_id\torf_id\tscore\tpval\tqval\textra_json" > ${meta.id}.orf_calls.tsv
  echo -e "chr1\t120\t210\torf|TX1|ribotricer\t0\t+\t120\t210\t0,0,0\t1\t90,\t0," > ${meta.id}.orf_calls.bed12
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      parser: stub
  END_VERSIONS
  """
}
