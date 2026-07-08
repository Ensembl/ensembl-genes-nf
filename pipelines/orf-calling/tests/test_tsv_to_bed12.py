import io
import csv
from pathlib import Path


def run_convert(tmp_path: Path, rows):
    tsv = tmp_path / "in.tsv"
    bed = tmp_path / "out.bed12"
    with tsv.open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample_id','tool','chrom','start','end','strand','frame','transcript_id','orf_id','score','pval','qval','extra_json'])
        for r in rows:
            w.writerow(r)
    # Import converter directly
    from importlib.util import spec_from_file_location, module_from_spec
    mod_path = Path(__file__).resolve().parents[1] / 'bin' / 'tsv_to_bed12.py'
    spec = spec_from_file_location('tsv_to_bed12', mod_path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore
    assert mod.convert(tsv, bed) >= 1
    return bed.read_text().strip().splitlines()


def test_bed12_genome_space(tmp_path):
    rows = [[
        'S1','ribotaper','chr1','100','160','+','','','ORF1','9','','','{}'
    ]]
    lines = run_convert(tmp_path, rows)
    assert len(lines) == 1
    f = lines[0].split('\t')
    assert f[0] == 'chr1' and f[1] == '100' and f[2] == '160' and f[5] == '+' and f[9] == '1'


def test_bed12_transcript_space(tmp_path):
    # No chrom -> falls back to transcript_id
    rows = [[
        'S1','ribocode','', '50','120','+','','TX1','ORF2','', '','', '{}'
    ]]
    lines = run_convert(tmp_path, rows)
    f = lines[0].split('\t')
    assert f[0] == 'TX1' and f[1] == '50' and f[2] == '120'

