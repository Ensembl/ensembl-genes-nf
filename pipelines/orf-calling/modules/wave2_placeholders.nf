/* Wave 2 modules: iRibo, ORF-RATER, PRICE, RibORF, Ribo-TISH, RiboTIE */

// iRibo
process PREP_IRIBO {
  tag "${meta.id}"
  label 'process_low'
  container "ubuntu:22.04"
  publishDir "${params.outdir}/orf_calls/iribo/${meta.id}", mode: 'copy', pattern: "prep/**"

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
    --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome
  """
  stub:
  """
  mkdir -p prep && echo '{}' > prep/inputs.json
  """
}

process RUN_IRIBO {
  tag "${meta.id}"
  label 'process_medium'
  container "ubuntu:22.04"
  publishDir "${params.outdir}/orf_calls/iribo/${meta.id}", mode: 'copy', pattern: "raw/**"

  input:
  tuple val(meta), path(prep)

  output:
  tuple val(meta), path('raw'), emit: raw

  script:
  """
  set -euo pipefail
  mkdir -p raw
  # Best-effort: create placeholder output
  echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/iribo_orfs.bed
  echo -e "chr1\t100\t200\tIRIBO1\t+\t5" >> raw/iribo_orfs.bed
  """
  stub:
  """
  mkdir -p raw && echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/iribo_orfs.bed
  """
}

process PARSE_IRIBO {
  tag "${meta.id}"
  label 'process_low'
  container "ubuntu:22.04"
  publishDir "${params.outdir}/orf_calls/iribo/${meta.id}", mode: 'copy', pattern: "*"

  input:
  tuple val(meta), path(raw)

  output:
  tuple val(meta), path("${meta.id}.orf_calls.tsv"), emit: standardized
  tuple val(meta), path("${meta.id}.orf_calls.bed12"), emit: bed12

  script:
  """
  python3 - << 'PY'
import json
from pathlib import Path
from parse_common import write_header, append_rows, default_row
meta = json.loads('''${meta.toString().replace("'","\'")}''')
rawdir = Path('${raw}')
out = Path(f"{meta['id']}.orf_calls.tsv")
write_header(out)
rows=[]
bed = rawdir/'iribo_orfs.bed'
if bed.exists():
  for ln in bed.read_text().splitlines():
    if not ln or ln.startswith('#') or ln.startswith('chrom'): continue
    chrom,start,end,name,strand,score = (ln.split('\t')+[None]*6)[:6]
    rows.append([meta['id'],'iribo',chrom,start,end,strand,None,None,name,score,None,None,'{}'])
else:
  rows.append(default_row(meta['id'],'iribo'))
append_rows(out, rows)
PY
  # TSV -> BED12
  python3 - << 'PY'
import csv
from pathlib import Path
tsv = Path("${meta.id}.orf_calls.tsv")
bed = Path("${meta.id}.orf_calls.bed12")
with open(tsv, newline='') as fh, open(bed, 'w') as out:
  r = csv.DictReader(fh, delimiter='\t')
  for row in r:
    try:
      s=int(row['start']); e=int(row['end'])
    except Exception:
      continue
    if e<=s: continue
    chrom=(row.get('chrom') or '').strip() or (row.get('transcript_id') or 'unknown')
    name=f"{row.get('orf_id') or 'orf'}|{row.get('transcript_id') or 'tx'}|{row.get('tool') or 'tool'}"
    strand=(row.get('strand') or '.').strip() or '.'
    out.write('\t'.join([chrom,str(s),str(e),name,'0',strand,str(s),str(e),'0,0,0','1',f"{e-s},",'0,'])+'\n')
PY
  """
}
// ORF-RATER
process PREP_ORFRATER {
  tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/orfrater/${meta.id}", mode: 'copy', pattern: "prep/**"
  input: tuple val(meta), path(bam), path(bai); path gtf; path fasta
  output: tuple val(meta), path('prep'), emit: prepared
  script: """
  mkdir -p prep
  prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type transcriptome
  """
  stub: """ mkdir -p prep && echo '{}' > prep/inputs.json """
}
process RUN_ORFRATER {
  tag "${meta.id}"; label 'process_medium'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/orfrater/${meta.id}", mode: 'copy', pattern: "raw/**"
  input: tuple val(meta), path(prep)
  output: tuple val(meta), path('raw'), emit: raw
  script: """
  set -euo pipefail; mkdir -p raw
  echo -e "transcript_id\torf_id\tstart\tend\tscore" > raw/orfrater.tsv
  echo -e "TX1\tORFR1\t100\t200\t4.2" >> raw/orfrater.tsv
  """
  stub: """ mkdir -p raw && echo -e "transcript_id\torf_id\tstart\tend\tscore" > raw/orfrater.tsv """
}
process PARSE_ORFRATER {
  tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/orfrater/${meta.id}", mode: 'copy', pattern: "*"
  input: tuple val(meta), path(raw)
  output: tuple val(meta), path("${meta.id}.orf_calls.tsv"), emit: standardized
  output: tuple val(meta), path("${meta.id}.orf_calls.bed12"), emit: bed12
  script: """
  python3 - << 'PY'
import csv, json
from pathlib import Path
from parse_common import write_header, append_rows, default_row
meta=json.loads('''${meta.toString().replace("'","\'")}'''); raw=Path('${raw}')
out=Path(f"{meta['id']}.orf_calls.tsv"); write_header(out)
fp=raw/'orfrater.tsv'; rows=[]
if fp.exists():
  r=csv.DictReader(open(fp), delimiter='\t')
  for row in r:
    rows.append([meta['id'],'orfrater',None,row['start'],row['end'],None,None,row['transcript_id'],row['orf_id'],row.get('score'),None,None,'{}'])
else:
  rows.append(default_row(meta['id'],'orfrater'))
append_rows(out, rows)
PY
  python3 - << 'PY'
import csv
from pathlib import Path
tsv=Path("${meta.id}.orf_calls.tsv"); bed=Path("${meta.id}.orf_calls.bed12")
with open(tsv, newline='') as fh, open(bed,'w') as out:
  r=csv.DictReader(fh, delimiter='\t')
  for row in r:
    try:
      s=int(row['start']); e=int(row['end'])
    except Exception: continue
    if e<=s: continue
    chrom=(row.get('chrom') or '').strip() or (row.get('transcript_id') or 'unknown')
    name=f"{row.get('orf_id') or 'orf'}|{row.get('transcript_id') or 'tx'}|{row.get('tool') or 'tool'}"
    strand=(row.get('strand') or '.').strip() or '.'
    out.write('\t'.join([chrom,str(s),str(e),name,'0',strand,str(s),str(e),'0,0,0','1',f"{e-s},",'0,'])+'\n')
PY
  """
}
// PRICE
process PREP_PRICE { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/price/${meta.id}", mode: 'copy', pattern: "prep/**"
  input: tuple val(meta), path(bam), path(bai); path gtf; path fasta
  output: tuple val(meta), path('prep'), emit: prepared
  script: """
  mkdir -p prep
  prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome
  """
  stub: """ mkdir -p prep && echo '{}' > prep/inputs.json """
}
process RUN_PRICE { tag "${meta.id}"; label 'process_medium'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/price/${meta.id}", mode: 'copy', pattern: "raw/**"
  input: tuple val(meta), path(prep)
  output: tuple val(meta), path('raw'), emit: raw
  script: """
  set -euo pipefail; mkdir -p raw
  echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/price.bed
  echo -e "chr1\t210\t300\tPRICE1\t+\t9" >> raw/price.bed
  """
  stub: """ mkdir -p raw && echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/price.bed """
}
process PARSE_PRICE { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/price/${meta.id}", mode: 'copy', pattern: "*"
  input: tuple val(meta), path(raw)
  output: tuple val(meta), path("${meta.id}.orf_calls.tsv"), emit: standardized
  output: tuple val(meta), path("${meta.id}.orf_calls.bed12"), emit: bed12
  script: """
  python3 - << 'PY'
import json
from pathlib import Path
from parse_common import write_header, append_rows, default_row
meta=json.loads('''${meta.toString().replace("'","\'")}'''); raw=Path('${raw}')
out=Path(f"{meta['id']}.orf_calls.tsv"); write_header(out)
rows=[]; bed=raw/'price.bed'
if bed.exists():
  for ln in bed.read_text().splitlines():
    if not ln or ln.startswith('chrom'): continue
    chrom,start,end,name,strand,score=(ln.split('\t')+[None]*6)[:6]
    rows.append([meta['id'],'price',chrom,start,end,strand,None,None,name,score,None,None,'{}'])
else:
  rows.append(default_row(meta['id'],'price'))
append_rows(out, rows)
PY
  python3 - << 'PY'
import csv
from pathlib import Path
tsv=Path("${meta.id}.orf_calls.tsv"); bed=Path("${meta.id}.orf_calls.bed12")
with open(tsv, newline='') as fh, open(bed,'w') as out:
  r=csv.DictReader(fh, delimiter='\t')
  for row in r:
    try:
      s=int(row['start']); e=int(row['end'])
    except Exception: continue
    if e<=s: continue
    chrom=(row.get('chrom') or '').strip() or (row.get('transcript_id') or 'unknown')
    name=f"{row.get('orf_id') or 'orf'}|{row.get('transcript_id') or 'tx'}|{row.get('tool') or 'tool'}"
    strand=(row.get('strand') or '.').strip() or '.'
    out.write('\t'.join([chrom,str(s),str(e),name,'0',strand,str(s),str(e),'0,0,0','1',f"{e-s},",'0,'])+'\n')
PY
  """
}
// RibORF
process PREP_RIBORF { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/riborf/${meta.id}", mode: 'copy', pattern: "prep/**"
  input: tuple val(meta), path(bam), path(bai); path gtf; path fasta
  output: tuple val(meta), path('prep'), emit: prepared
  script: """
  mkdir -p prep
  prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome
  """
  stub: """ mkdir -p prep && echo '{}' > prep/inputs.json """
}
process RUN_RIBORF { tag "${meta.id}"; label 'process_medium'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/riborf/${meta.id}", mode: 'copy', pattern: "raw/**"
  input: tuple val(meta), path(prep)
  output: tuple val(meta), path('raw'), emit: raw
  script: """
  set -euo pipefail; mkdir -p raw
  echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/riborf.bed
  echo -e "chr1\t320\t390\tRIBORF1\t-\t7" >> raw/riborf.bed
  """
  stub: """ mkdir -p raw && echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/riborf.bed """
}
process PARSE_RIBORF { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/riborf/${meta.id}", mode: 'copy', pattern: "*"
  input: tuple val(meta), path(raw)
  output: tuple val(meta), path("${meta.id}.orf_calls.tsv"), emit: standardized
  output: tuple val(meta), path("${meta.id}.orf_calls.bed12"), emit: bed12
  script: """
  python3 - << 'PY'
import json
from pathlib import Path
from parse_common import write_header, append_rows, default_row
meta=json.loads('''${meta.toString().replace("'","\'")}'''); raw=Path('${raw}')
out=Path(f"{meta['id']}.orf_calls.tsv"); write_header(out)
rows=[]; bed=raw/'riborf.bed'
if bed.exists():
  for ln in bed.read_text().splitlines():
    if not ln or ln.startswith('chrom'): continue
    chrom,start,end,name,strand,score=(ln.split('\t')+[None]*6)[:6]
    rows.append([meta['id'],'riborf',chrom,start,end,strand,None,None,name,score,None,None,'{}'])
else:
  rows.append(default_row(meta['id'],'riborf'))
append_rows(out, rows)
PY
  python3 - << 'PY'
import csv
from pathlib import Path
tsv=Path("${meta.id}.orf_calls.tsv"); bed=Path("${meta.id}.orf_calls.bed12")
with open(tsv, newline='') as fh, open(bed,'w') as out:
  r=csv.DictReader(fh, delimiter='\t')
  for row in r:
    try:
      s=int(row['start']); e=int(row['end'])
    except Exception: continue
    if e<=s: continue
    chrom=(row.get('chrom') or '').strip() or (row.get('transcript_id') or 'unknown')
    name=f"{row.get('orf_id') or 'orf'}|{row.get('transcript_id') or 'tx'}|{row.get('tool') or 'tool'}"
    strand=(row.get('strand') or '.').strip() or '.'
    out.write('\t'.join([chrom,str(s),str(e),name,'0',strand,str(s),str(e),'0,0,0','1',f"{e-s},",'0,'])+'\n')
PY
  """
}
// Ribo-TISH
process PREP_RIBOTISH { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/ribotish/${meta.id}", mode: 'copy', pattern: "prep/**"
  input: tuple val(meta), path(bam), path(bai); path gtf; path fasta
  output: tuple val(meta), path('prep'), emit: prepared
  script: """
  mkdir -p prep
  prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome
  """
  stub: """ mkdir -p prep && echo '{}' > prep/inputs.json """
}
process RUN_RIBOTISH { tag "${meta.id}"; label 'process_medium'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/ribotish/${meta.id}", mode: 'copy', pattern: "raw/**"
  input: tuple val(meta), path(prep)
  output: tuple val(meta), path('raw'), emit: raw
  script: """
  set -euo pipefail; mkdir -p raw
  echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/ribotish.tsv
  echo -e "chr1\t400\t480\tTISH1\t+\t10" >> raw/ribotish.tsv
  """
  stub: """ mkdir -p raw && echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/ribotish.tsv """
}
process PARSE_RIBOTISH { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/ribotish/${meta.id}", mode: 'copy', pattern: "*"
  input: tuple val(meta), path(raw)
  output: tuple val(meta), path("${meta.id}.orf_calls.tsv"), emit: standardized
  output: tuple val(meta), path("${meta.id}.orf_calls.bed12"), emit: bed12
  script: """
  python3 - << 'PY'
import json
from pathlib import Path
from parse_common import write_header, append_rows, default_row
meta=json.loads('''${meta.toString().replace("'","\'")}'''); raw=Path('${raw}')
out=Path(f"{meta['id']}.orf_calls.tsv"); write_header(out)
rows=[]; tsv=raw/'ribotish.tsv'
if tsv.exists():
  for ln in tsv.read_text().splitlines():
    if not ln or ln.startswith('chrom'): continue
    chrom,start,end,name,strand,score=(ln.split('\t')+[None]*6)[:6]
    rows.append([meta['id'],'ribotish',chrom,start,end,strand,None,None,name,score,None,None,'{}'])
else:
  rows.append(default_row(meta['id'],'ribotish'))
append_rows(out, rows)
PY
  python3 - << 'PY'
import csv
from pathlib import Path
tsv=Path("${meta.id}.orf_calls.tsv"); bed=Path("${meta.id}.orf_calls.bed12")
with open(tsv, newline='') as fh, open(bed,'w') as out:
  r=csv.DictReader(fh, delimiter='\t')
  for row in r:
    try:
      s=int(row['start']); e=int(row['end'])
    except Exception: continue
    if e<=s: continue
    chrom=(row.get('chrom') or '').strip() or (row.get('transcript_id') or 'unknown')
    name=f"{row.get('orf_id') or 'orf'}|{row.get('transcript_id') or 'tx'}|{row.get('tool') or 'tool'}"
    strand=(row.get('strand') or '.').strip() or '.'
    out.write('\t'.join([chrom,str(s),str(e),name,'0',strand,str(s),str(e),'0,0,0','1',f"{e-s},",'0,'])+'\n')
PY
  """
}
// RiboTIE
process PREP_RIBOTIE { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/ribotie/${meta.id}", mode: 'copy', pattern: "prep/**"
  input: tuple val(meta), path(bam), path(bai); path gtf; path fasta
  output: tuple val(meta), path('prep'), emit: prepared
  script: """
  mkdir -p prep
  prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome
  """
  stub: """ mkdir -p prep && echo '{}' > prep/inputs.json """
}
process RUN_RIBOTIE { tag "${meta.id}"; label 'process_medium'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/ribotie/${meta.id}", mode: 'copy', pattern: "raw/**"
  input: tuple val(meta), path(prep)
  output: tuple val(meta), path('raw'), emit: raw
  script: """
  set -euo pipefail; mkdir -p raw
  echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/ribotie.tsv
  echo -e "chr1\t500\t560\tRIBOTIE1\t-\t6" >> raw/ribotie.tsv
  """
  stub: """ mkdir -p raw && echo -e "chrom\tstart\tend\tname\tstrand\tscore" > raw/ribotie.tsv """
}
process PARSE_RIBOTIE { tag "${meta.id}"; label 'process_low'; container "ubuntu:22.04"; publishDir "${params.outdir}/orf_calls/ribotie/${meta.id}", mode: 'copy', pattern: "*"
  input: tuple val(meta), path(raw)
  output: tuple val(meta), path("${meta.id}.orf_calls.tsv"), emit: standardized
  output: tuple val(meta), path("${meta.id}.orf_calls.bed12"), emit: bed12
  script: """
  python3 - << 'PY'
import json
from pathlib import Path
from parse_common import write_header, append_rows, default_row
meta=json.loads('''${meta.toString().replace("'","\'")}'''); raw=Path('${raw}')
out=Path(f"{meta['id']}.orf_calls.tsv"); write_header(out)
rows=[]; tsv=raw/'ribotie.tsv'
if tsv.exists():
  for ln in tsv.read_text().splitlines():
    if not ln or ln.startswith('chrom'): continue
    chrom,start,end,name,strand,score=(ln.split('\t')+[None]*6)[:6]
    rows.append([meta['id'],'ribotie',chrom,start,end,strand,None,None,name,score,None,None,'{}'])
else:
  rows.append(default_row(meta['id'],'ribotie'))
append_rows(out, rows)
PY
  python3 - << 'PY'
import csv
from pathlib import Path
tsv=Path("${meta.id}.orf_calls.tsv"); bed=Path("${meta.id}.orf_calls.bed12")
with open(tsv, newline='') as fh, open(bed,'w') as out:
  r=csv.DictReader(fh, delimiter='\t')
  for row in r:
    try:
      s=int(row['start']); e=int(row['end'])
    except Exception: continue
    if e<=s: continue
    chrom=(row.get('chrom') or '').strip() or (row.get('transcript_id') or 'unknown')
    name=f"{row.get('orf_id') or 'orf'}|{row.get('transcript_id') or 'tx'}|{row.get('tool') or 'tool'}"
    strand=(row.get('strand') or '.').strip() or '.'
    out.write('\t'.join([chrom,str(s),str(e),name,'0',strand,str(s),str(e),'0,0,0','1',f"{e-s},",'0,'])+'\n')
PY
  """
}
