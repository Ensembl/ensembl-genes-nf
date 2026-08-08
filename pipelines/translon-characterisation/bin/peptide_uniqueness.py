#!/usr/bin/env python3
"""Classify exact/near-exact GENCODE proteome overlap from MMseqs2 results."""
import argparse, json, subprocess, tempfile
from pathlib import Path

def read(path):
    with open(path) as h: yield from (json.loads(x) for x in h if x.strip())

def main():
    p=argparse.ArgumentParser(); p.add_argument('--peptides',required=True); p.add_argument('--proteome',required=True); p.add_argument('--output',required=True); p.add_argument('--mmseqs',default='mmseqs'); a=p.parse_args()
    products=list(read(a.peptides)); valid=[x for x in products if x.get('sequence')]
    hits={}
    if valid:
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); query=root/'query.fa'; query.write_text(''.join(f">{x['peptide_id']}\n{x['sequence']}\n" for x in valid))
            qdb, pdb, result, work = (root/'qdb', root/'pdb', root/'result', root/'work')
            for cmd in ([a.mmseqs,'createdb',str(query),str(qdb)], [a.mmseqs,'createdb',a.proteome,str(pdb)], [a.mmseqs,'search',str(qdb),str(pdb),str(result),str(work),'--min-seq-id','0.95','-c','1.0'], [a.mmseqs,'convertalis',str(qdb),str(pdb),str(result),str(root/'hits.tsv'),'--format-output','query,target,pident,qcov,tcov']): subprocess.run(cmd,check=True)
            for line in (root/'hits.tsv').read_text().splitlines():
                q,t,pid,qcov,tcov=line.split('\t'); hits.setdefault(q,[]).append((t,float(pid),float(qcov),float(tcov)))
    with open(a.output,'w') as out:
        for x in products:
            if not x.get('sequence'):
                x['peptide_uniqueness']={'state':'unattributable','reason':'no trusted translated peptide sequence'}
            elif not hits.get(x['peptide_id']):
                x['peptide_uniqueness']={'state':'positive','class':'novel','reason':'no >=95% identity full-query GENCODE match'}
            else:
                target,pid,qcov,tcov=sorted(hits[x['peptide_id']],key=lambda h:(h[1],h[2],h[3]),reverse=True)[0]
                cls='identical' if tcov >= .999 else 'substring'
                x['peptide_uniqueness']={'state':'positive','class':cls,'parent_protein_id':target,'identity':pid,'query_coverage':qcov,'target_coverage':tcov}
            out.write(json.dumps(x,sort_keys=True)+'\n')
if __name__=='__main__': main()
