#!/usr/bin/env python3
import argparse, json
def main():
    p=argparse.ArgumentParser(); p.add_argument('--instances',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    with open(a.output,'w') as out:
        for x in (json.loads(line) for line in open(a.instances) if line.strip()):
            parts=[]
            for exon in x.get('exon_blocks',[]):
                start=max(x['start'],exon['start']); end=min(x['end'],exon['end'])
                if start < end: parts.append(f"{x['chrom']}:{start+1}-{end}")
            if parts: out.write(f"{'+'.join(parts)}\t{x['strand']}\t{x['instance_id']}\n")
if __name__=='__main__': main()
