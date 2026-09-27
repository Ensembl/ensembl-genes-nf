#!/usr/bin/env python3
"""Emit mafExtract BED6 blocks from exact exon intersections."""
import argparse, json, re
from pathlib import Path

def main():
    p=argparse.ArgumentParser(); p.add_argument('--instances',required=True); p.add_argument('--outdir',required=True); p.add_argument('--manifest',required=True); a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); manifest=[]
    for record in (json.loads(line) for line in open(a.instances) if line.strip()):
        blocks=record.get('exon_blocks')
        if not blocks: continue
        safe=re.sub(r'[^A-Za-z0-9_.-]+','_',record['instance_id']); paths=[]
        for n,block in enumerate(blocks,1):
            left=max(record['start'],block['start']); right=min(record['end'],block['end'])
            if left >= right: continue
            path=out/f'{safe}.exon{n}.bed'
            path.write_text(f"{record['chrom']}\t{left}\t{right}\t{safe}_exon{n}\t0\t{record['strand']}\n")
            paths.append(str(path))
        if paths: manifest.append({'instance_id':record['instance_id'],'chrom':record['chrom'],'strand':record['strand'],'beds':paths})
    with open(a.manifest,'w') as h:
        for row in manifest: h.write(json.dumps(row,sort_keys=True)+'\n')
if __name__=='__main__': main()
