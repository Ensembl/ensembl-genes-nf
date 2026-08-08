#!/usr/bin/env python3
import argparse, json
def rows(path):
    with open(path) as h: yield from (json.loads(x) for x in h if x.strip())
def main():
    p=argparse.ArgumentParser(); p.add_argument('--instances',required=True); p.add_argument('--peptides',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    by_id={x['peptide_id']: x.get('peptide_uniqueness', {'state':'unattributable','reason':'missing peptide result'}) for x in rows(a.peptides)}
    with open(a.output,'w') as out:
        for x in rows(a.instances):
            key=f"{x['interval_id']}|{x.get('translation',{}).get('frame',x.get('frame','.'))}"
            x['axes']=x.get('axes',{}); x['axes']['peptide_uniqueness']=by_id.get(key,{'state':'unattributable','reason':'no peptide product for instance'})
            out.write(json.dumps(x,sort_keys=True)+'\n')
if __name__=='__main__': main()
