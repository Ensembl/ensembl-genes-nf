#!/usr/bin/env python3
import argparse, json
def rows(path):
    with open(path) as h: yield from (json.loads(x) for x in h if x.strip())
def main():
    p=argparse.ArgumentParser(); p.add_argument('--instances',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    with open(a.output,'w') as out:
        for x in rows(a.instances):
            cls=x['class']; regulatory=cls in {'uORF','uoORF','dORF','doORF'}
            state='positive' if regulatory else 'null'
            x['axes']=x.get('axes',{}); x['axes']['regulatory_geometry']={'state':state,'reason':'CDS-relative topology supports a regulatory geometry' if regulatory else 'no CDS-relative regulatory geometry','class':cls,'effect_direction':'unresolved' if regulatory else None,'attribution_target':'translation_act'}
            out.write(json.dumps(x,sort_keys=True)+'\n')
if __name__=='__main__': main()
