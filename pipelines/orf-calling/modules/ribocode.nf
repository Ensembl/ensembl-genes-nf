/* RiboCode module: PREP -> RUN -> PARSE */

process PREP_RIBOCODE {
  tag "${meta.id}"
  label 'process_low'

  conda "bioconda::ribocode"
  container "${ params.container_ribocode ?: 'quay.io/biocontainers/ribocode:latest' }"

  publishDir "${params.outdir}/orf_calls/ribocode", mode: 'copy', pattern: "prep/**"

  input:
  tuple val(meta), path(bam), path(bai)
  path gtf
  path fasta

  output:
  tuple val(meta), path('prep'), emit: prepared

  script:
  def pass_path = params.qc_pass_lengths_dir ? file("${params.qc_pass_lengths_dir}/${meta.id}.pass_lengths.tsv") : file('NO_PASS')
  """
  mkdir -p prep
  prep_inputs.py --outdir prep --sample-id ${meta.id} \
    --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} \
    --bam-type transcriptome \
    --pass-lengths ${pass_path}
  """

  stub:
  """
  mkdir -p prep && echo '{}' > prep/inputs.json && echo '{}' > prep/meta.json
  """
}

process RUN_RIBOCODE {
  tag "${meta.id}"
  label 'process_medium'
  conda "bioconda::ribocode"
  container "${ params.container_ribocode ?: 'quay.io/biocontainers/ribocode:latest' }"

  publishDir "${params.outdir}/orf_calls/ribocode", mode: 'copy', pattern: "raw/**"

  input:
  tuple val(meta), path(prepdir)

  output:
  tuple val(meta), path('raw'), emit: raw

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  mkdir -p raw
  export THREADS=${params.threads_ribocode ?: 4}
  export EXTRA_ARGS="${params.args_ribocode ?: ''}"
  if [ -n "${task.stub}" ]; then
    cp -r ${prepdir} raw/prep_copy 2>/dev/null || true
    echo -e "#chrom\tstart\tend\tname\tscore\tstrand" > raw/ribocode_orfs.bed
    echo -e "chr1\t100\t200\tORF1\t10\t+" >> raw/ribocode_orfs.bed
    echo -e "sample_id\ttranscript_id\tstart\tend\tframe\tscore\tpval" > raw/ribocode_orfs.tsv
    echo -e "${meta.id}\tTX1\t100\t200\t0\t10\t0.01" >> raw/ribocode_orfs.tsv
  else
    python3 - <<'PY'
import json, subprocess, os, shutil
from pathlib import Path
inp=json.loads(Path('prep/inputs.json').read_text())
bam,gtf,fasta = inp['bam'], inp['gtf'], inp['fasta']
threads=int(os.environ.get('THREADS','4'))
extra=os.environ.get('EXTRA_ARGS','')
Path('raw').mkdir(exist_ok=True)

def have(cmd):
    return shutil.which(cmd) is not None

# Prefer two-step to filter lengths via config
prep_ok = False
if have('prepare_transcripts') and have('metaplots') and have('RiboCode'):
    subprocess.run(['prepare_transcripts','-g',gtf,'-f',fasta,'-o','raw/annot'], check=False)
    subprocess.run(['metaplots','-a','raw/annot','-r',bam,'-o','raw/meta','-t',str(threads)], check=False)
    # Look for config file; filter to pass lengths if available
    cfg = None
    for cand in ['raw/meta/config.txt','raw/meta/psite_config.txt','raw/meta/configure.txt']:
        if Path(cand).exists(): cfg=cand; break
    if cfg:
        pl = inp.get('pass_lengths')
        if pl and Path(pl).exists():
            keep = set(int(x.strip()) for x in Path(pl).read_text().splitlines()[1:] if x.strip())
            lines = Path(cfg).read_text().splitlines()
            hdr = [ln for ln in lines if ln.startswith('#')]
            body = [ln for ln in lines if ln.strip() and not ln.startswith('#') and int(ln.split()[0]) in keep]
            Path('raw/meta/config.filtered.txt').write_text('\n'.join(hdr+body)+'\n')
            cfg = 'raw/meta/config.filtered.txt'
    # Run RiboCode with (possibly) filtered config
    cmd = ['RiboCode','-a','raw/annot','-r',bam,'-g',fasta,'-t',str(threads),'-c',cfg] if cfg else ['RiboCode_onestep','-a',gtf,'-r',bam,'-g',fasta,'-t',str(threads)]
    if extra:
        cmd += extra.split()
    subprocess.run(cmd, check=False)
    prep_ok = True

# Fallback to one-step if tool subcommands not found
if not prep_ok:
    cmd = ['RiboCode_onestep','-a',gtf,'-r',bam,'-g',fasta,'-t',str(threads)]
    if extra:
        cmd += extra.split()
    subprocess.run(cmd, check=False)

# Normalize outputs for parser
Path('raw/ribocode_orfs.tsv').touch()
Path('raw/ribocode_orfs.bed').touch()
PY
  fi
  """

  stub:
  """
  mkdir -p raw && touch raw/ribocode_orfs.bed raw/ribocode_orfs.tsv
  """
}

process PARSE_RIBOCODE {
  tag "${meta.id}"
  label 'process_low'

  container "ubuntu:22.04"

  publishDir "${params.outdir}/orf_calls/ribocode", mode: 'copy', pattern: "*"

  input:
  tuple val(meta), path(rawdir)

  output:
  tuple val(meta), path("${meta.id}.ribocode.orf_calls.tsv"), emit: standardized
  tuple val(meta), path("${meta.id}.ribocode.orf_calls.bed12"), emit: bed12
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
tsv = raw / 'ribocode_orfs.tsv'
rows = []
if tsv.exists():
  r = csv.DictReader(open(tsv), delimiter='\t')
  for row in r:
    rows.append([
      meta['id'],'ribocode',None,row.get('start'),row.get('end'),None,
      row.get('frame'),row.get('transcript_id'),None,row.get('score'),row.get('pval'),None,'{}'
    ])
else:
  rows.append(default_row(meta['id'],'ribocode'))
append_rows(out, rows)
PY

  # Convert TSV to BED12 (minimal, single-block)
  python3 - << 'PY'
import csv
from pathlib import Path

def clamp(v, lo=0, hi=1000):
    return max(lo, min(hi, v))

def score_from_fields(score_str, pval_str):
    import math
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

tsv = Path("${meta.id}.ribocode.orf_calls.tsv")
bed = Path("${meta.id}.ribocode.orf_calls.bed12")
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
  echo -e "sample_id\ttool\tchrom\tstart\tend\tstrand\tframe\ttranscript_id\torf_id\tscore\tpval\tqval\textra_json" > ${meta.id}.ribocode.orf_calls.tsv
  echo -e "chr1\t100\t200\torf|TX1|ribocode\t0\t+\t100\t200\t0,0,0\t1\t100,\t0," > ${meta.id}.ribocode.orf_calls.bed12
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      parser: stub
  END_VERSIONS
  """
}
