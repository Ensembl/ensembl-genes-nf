#!/usr/bin/env python3
"""Wrap a raw PhyloCSF result with alignment-depth provenance.

Interpretation against length/GC/class/depth-matched nulls remains a required
upstream contract; this script will not call a raw score positive by itself.
"""
import argparse, json, re
from pathlib import Path


def alignment_instance_id(path):
    """Read the identity written into the downloaded MSA reference header.

    ORBL_tools file names are an implementation detail, so fan-back must not
    infer identity from them.  The download normaliser writes this header into
    the FASTA itself before the file enters a per-alignment process.
    """
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                fields = line[1:].strip().split()
                for field in fields[1:]:
                    if field.startswith('instance_id='):
                        return field.split('=', 1)[1]
                break
    return None


def main():
    p=argparse.ArgumentParser(); p.add_argument('--raw',required=True); p.add_argument('--output',required=True); p.add_argument('--matched-null',default=None); p.add_argument('--alignment'); p.add_argument('--instance-id'); p.add_argument('--identities'); a=p.parse_args()
    if a.identities:
        identities = {Path(x['alignment']).name: x['instance_id'] for x in (json.loads(line) for line in open(a.identities) if line.strip())}
        null = json.load(open(a.matched_null)) if a.matched_null else None
        results = []
        for line in open(a.raw):
            fields = line.rstrip().split('\t')
            if len(fields) < 2:
                continue
            name, status, value = Path(fields[0]).name, fields[1], fields[-1]
            iid = identities.get(name)
            match = re.fullmatch(r'-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', value)
            common = {'axis': 'phylocsf'}
            if not iid:
                result = {**common, 'state': 'unattributable', 'reason': 'PhyloCSF alignment has no identity mapping', 'alignment': name}
            elif status == 'abort' or not match:
                result = {**common, 'state': 'uninformative', 'reason': 'PhyloCSF did not emit a score', 'status': status}
            elif null is None:
                result = {**common, 'state': 'unattributable', 'reason': 'raw score lacks a matched null distribution', 'raw_score': float(value)}
            else:
                score, threshold = float(value), float(null['positive_threshold'])
                result = {**common, 'state': 'positive' if score >= threshold else 'null', 'reason': 'matched-null comparison', 'raw_score': score, 'positive_threshold': threshold, 'null_provenance': null.get('provenance'), 'null_match': {key: null.get(key) for key in ('length_bin', 'gc_bin', 'class', 'depth_bin', 'reference', 'model')}}
            results.append({'instance_id': iid, 'axis': {'phylocsf': result}} if iid else {'axis': {'phylocsf': result}})
        Path(a.output).write_text(''.join(json.dumps(x, sort_keys=True) + '\n' for x in results))
        return
    text=open(a.raw).read(); match=re.search(r'(-?\d+(?:\.\d+)?)',text)
    instance_id = a.instance_id or (alignment_instance_id(a.alignment) if a.alignment else None)
    common = {'axis': 'phylocsf'}
    if instance_id:
        common['instance_id'] = instance_id
    if not match:
        result={**common, 'state':'uninformative','reason':'PhyloCSF emitted no parseable score'}
    elif not a.matched_null:
        result={**common, 'state':'unattributable','reason':'raw score lacks a matched null distribution','raw_score':float(match.group(1))}
    else:
        null=json.load(open(a.matched_null)); score=float(match.group(1)); threshold=float(null['positive_threshold'])
        result={**common, 'state':'positive' if score>=threshold else 'null','reason':'matched-null comparison','raw_score':score,'positive_threshold':threshold,'null_provenance':null.get('provenance'), 'null_match': {key: null.get(key) for key in ('length_bin', 'gc_bin', 'class', 'depth_bin', 'reference', 'model')}}
    if instance_id:
        json.dump({'instance_id': instance_id, 'axis': {'phylocsf': result}}, open(a.output, 'w'), sort_keys=True)
    else:
        # Do not permit a score with unknown attribution to fan back to every
        # transcript context.  The caller can retain this as diagnostic output.
        json.dump({'axis': {'phylocsf': {**result, 'state': 'unattributable', 'reason': 'alignment has no instance identity'}}}, open(a.output, 'w'), sort_keys=True)
if __name__=='__main__': main()
