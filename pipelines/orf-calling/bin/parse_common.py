#!/usr/bin/env python3
import argparse, csv, json, sys
from pathlib import Path

COLUMNS = [
    'sample_id','tool','chrom','start','end','strand','frame',
    'transcript_id','orf_id','score','pval','qval','extra_json'
]

def write_header(path):
    with open(path,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(COLUMNS)

def append_rows(path, rows):
    with open(path,'a',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerows(rows)

def default_row(sample_id, tool):
    return [sample_id, tool, None, None, None, None, None, None, None, None, None, None, '{}']

def load_meta(meta_path):
    try:
        return json.loads(Path(meta_path).read_text())
    except Exception:
        return {'id':'unknown'}

def safe_int(x):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return None

