#!/usr/bin/env python3
"""Map official ORBL_tools tabular output to the typed axis contract."""
import argparse, json
def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    with open(a.output,'w') as out:
        for line in open(a.input):
            if not line.strip() or line.startswith('#'): continue
            fields=line.rstrip('\n').split('\t')
            # input: intervals, strand, biotype, instance_id; output then v, q, start, stop, frame
            if len(fields) < 6: continue
            iid, v, q = fields[3], fields[4], fields[5]
            state = 'uninformative' if v == 'NA' else ('positive' if q not in ('NA','') and float(q) >= .9 else 'null' if q not in ('NA','') else 'unattributable')
            reason = 'ORBLq unavailable for this class/alignment set' if q == 'NA' else 'ORBL matched-null constraint score'
            record = {'instance_id': iid, 'axis': {'orbl': {'state': state, 'reason': reason, 'orblv': None if v == 'NA' else float(v), 'orblq': None if q == 'NA' else float(q), 'components': fields[6:]}}}
            out.write(json.dumps(record, sort_keys=True) + '\n')
if __name__=='__main__': main()
