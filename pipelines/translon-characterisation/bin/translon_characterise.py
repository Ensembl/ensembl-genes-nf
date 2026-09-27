#!/usr/bin/env python3
"""Pure-Python attribution and verdict layer for translon characterisation.

The workflow deliberately keeps biological decision logic here rather than in
Nextflow.  Inputs/outputs are tabular JSONL so each commodity-tool wrapper can
be independently replaced without changing adjudication semantics.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, re
from collections import defaultdict
from pathlib import Path

STATES = {'positive', 'null', 'uninformative', 'unattributable'}

def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)

def rows(path):
    with open_text(path) as handle:
        first = handle.read(1); handle.seek(0)
        if first == '{':
            yield from (json.loads(line) for line in handle if line.strip())
        else:
            yield from csv.DictReader(handle, delimiter='\t')

def write_rows(path, values):
    with open(path, 'w') as handle:
        for value in values:
            handle.write(json.dumps(value, sort_keys=True) + '\n')

def typed(state, reason, **extra):
    assert state in STATES
    return {'state': state, 'reason': reason, **extra}

def parse_intervals(path):
    """Read headered TSV/JSONL or standard BED3-BED12.

    BED name is the stable ORF identifier.  BED thickStart is not treated as an
    ORF start: these records already describe annotated ORFs, not display-only
    transcript models.  BED12 blocks are retained for provenance and checked.
    """
    with open_text(path) as handle:
        first = next((line for line in handle if line.strip() and not line.startswith('#')), '')
    fields = first.rstrip('\n').split('\t')
    headered = first.startswith('{') or ({'start', 'end'} <= set(fields))
    if headered:
        for n, row in enumerate(rows(path), 1):
            start, end = int(row['start']), int(row['end'])
            if start < 0 or end <= start:
                raise ValueError(f'invalid interval coordinates on record {n}: {start}-{end}')
            chrom = row.get('chrom') or row.get('seqid')
            if not chrom:
                raise ValueError(f'missing chromosome on record {n}')
            yield {'interval_id': row.get('interval_id') or row.get('name') or f'interval_{n}',
                   'chrom': chrom, 'start': start, 'end': end,
                   'strand': row.get('strand', '.'), 'frame': row.get('frame', '.')}
        return
    with open_text(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            bed = line.rstrip('\n').split('\t')
            if len(bed) < 3:
                raise ValueError(f'BED record {n} has fewer than 3 columns')
            chrom, start, end = bed[0], int(bed[1]), int(bed[2])
            if start < 0 or end <= start:
                raise ValueError(f'invalid BED coordinates on line {n}: {start}-{end}')
            strand = bed[5] if len(bed) >= 6 else '.'
            if strand not in ('+', '-', '.'):
                raise ValueError(f'invalid BED strand on line {n}: {strand}')
            item = {'interval_id': bed[3] if len(bed) >= 4 and bed[3] else f'interval_{n}',
                    'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
                    'frame': bed[12] if len(bed) >= 13 and bed[12] else '.'}
            if len(bed) >= 12:
                count = int(bed[9])
                sizes = [int(x) for x in bed[10].rstrip(',').split(',') if x]
                starts = [int(x) for x in bed[11].rstrip(',').split(',') if x]
                if count != len(sizes) or count != len(starts):
                    raise ValueError(f'BED12 block count mismatch on line {n}')
                blocks = [{'start': start + offset, 'end': start + offset + size}
                          for size, offset in zip(sizes, starts)]
                if any(b['start'] < start or b['end'] > end or b['end'] <= b['start'] for b in blocks):
                    raise ValueError(f'invalid BED12 blocks on line {n}')
                item['input_blocks'] = blocks
            yield item

def gff_features(path):
    genes, transcripts, exons, cds = {}, {}, defaultdict(list), defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9: continue
            attrs = {k: v for k, v in re.findall(r'([^=;]+)=([^;]+)', f[8])}
            ident = attrs.get('ID', '').split(':')[-1]
            parents = [x.split(':')[-1] for x in attrs.get('Parent', '').split(',') if x]
            item = {'chrom': f[0], 'type': f[2], 'start': int(f[3]) - 1, 'end': int(f[4]),
                    'strand': f[6], 'phase': f[7], 'attrs': attrs, 'id': ident, 'parents': parents}
            if f[2] in ('mRNA', 'transcript'):
                transcripts[ident] = item
            elif f[2] == 'exon':
                for p in parents: exons[p].append(item)
            elif f[2] == 'CDS':
                for p in parents: cds[p].append(item)
    return transcripts, exons, cds

def genome(path):
    result, name, sequence = {}, None, []
    for line in open_text(path):
        line = line.strip()
        if line.startswith('>'):
            if name: result[name] = ''.join(sequence).upper()
            name = line[1:].split()[0]; sequence = []
        elif line: sequence.append(line)
    if name: result[name] = ''.join(sequence).upper()
    return result

def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]

def spliced_orf_sequence(interval, exon_parts, seq):
    pieces = [(max(interval['start'], p['start']), min(interval['end'], p['end'])) for p in exon_parts]
    pieces = [(left, right) for left, right in pieces if left < right]
    pieces.sort(reverse=interval['strand'] == '-')
    dna = ''.join(seq[left:right] for left, right in pieces)
    return reverse_complement(dna) if interval['strand'] == '-' else dna

def overlap(a, b): return a['chrom'] == b['chrom'] and a['start'] < b['end'] and b['start'] < a['end']
def spans(parts, start, end):
    """True when every base of a genomic interval is in spliced exons."""
    covered = sorted((max(start,p['start']), min(end,p['end'])) for p in parts if max(start,p['start']) < min(end,p['end']))
    cursor = start
    for left, right in covered:
        if left > cursor: return False
        cursor = max(cursor, right)
    return cursor == end

def substrate(args):
    txs, exons, cds = gff_features(args.gff)
    sequences = genome(args.genome)
    instances, unhosted = [], []
    for interval in parse_intervals(args.intervals):
        compatible = []
        for tid, tx in txs.items():
            if not overlap(interval, tx) or interval['strand'] not in ('.', tx['strand']): continue
            if not spans(exons[tid], interval['start'], interval['end']): continue
            if interval['chrom'] not in sequences: continue
            candidate = {**interval, 'strand': tx['strand']}
            orf_dna = spliced_orf_sequence(candidate, exons[tid], sequences[interval['chrom']])
            # A non-AUG call can be admitted only by the trusted translation adapter later.
            strict_orf = len(orf_dna) >= 6 and len(orf_dna) % 3 == 0 and orf_dna[:3] == 'ATG' and orf_dna[-3:] in {'TAA','TAG','TGA'} and not any(orf_dna[i:i+3] in {'TAA','TAG','TGA'} for i in range(3, len(orf_dna)-3, 3))
            if not strict_orf: continue
            # This is a conservative compatibility gate. Start/stop/frame are reconciled later
            # against the trusted TranslonScorer frame because intervals carry coordinates only.
            overlaps_cds = any(overlap(interval, part) for part in cds[tid])
            class_ = 'intORF' if overlaps_cds else 'uORF'
            if cds[tid]:
                first = min(p['start'] for p in cds[tid]); last = max(p['end'] for p in cds[tid])
                if interval['end'] <= first: class_ = 'uORF'
                elif interval['start'] < first < interval['end']: class_ = 'uoORF'
                elif interval['start'] >= last: class_ = 'dORF'
            iid = f"{interval['interval_id']}|{tid}"
            instance = {**interval, 'instance_id': iid, 'transcript_id': tid,
                        'class': class_, 'overlaps_cds': overlaps_cds,
                        'exon_blocks': [{'start': p['start'], 'end': p['end']} for p in sorted(exons[tid], key=lambda p: p['start'])],
                        'orf_dna': orf_dna,
                        'admissibility': typed('positive', 'spliced exon-compatible transcript context')}
            compatible.append(instance)
        if compatible: instances.extend(compatible)
        else: unhosted.append({**interval, 'claim_hint': 'translated-no-compatible-transcript',
                               'admissibility': typed('unattributable', 'no compatible annotated transcript')})
    write_rows(args.instances, instances)
    write_rows(args.unhosted, unhosted)

def translation_join(args):
    verdicts = {r.get('interval_id') or r.get('id'): r for r in rows(args.verdicts)}
    out = []
    for instance in rows(args.instances):
        verdict = verdicts.get(instance['interval_id'], {})
        required = ('frame', 'mechanism_class', 'confidence', 'condition_state', 'peptide_sequence')
        missing = [key for key in required if verdict.get(key) in (None, '')]
        if missing:
            instance['translation'] = {'validation': 'failed', 'missing_fields': missing}
            instance['admissibility'] = typed('unattributable', 'trusted TranslonScorer record is incomplete', missing_fields=missing)
            out.append(instance)
            continue
        frame = verdict.get('frame', '.')
        agree = frame in ('.', '', instance.get('frame', '.'), verdict.get('transcript_frame', frame))
        instance['translation'] = {'frame': frame, 'mechanism': verdict.get('mechanism_class', 'unknown'),
                                   'confidence': verdict.get('confidence', 'unknown'), 'frame_agreement': agree,
                                   'condition_state': verdict.get('condition_state', 'condition-unresolved'),
                                   'biotype_with_frameshift': verdict.get('biotype_with_frameshift', instance['class'])}
        instance['peptide_sequence'] = verdict['peptide_sequence']
        if not agree: instance['admissibility'] = typed('unattributable', 'trusted frame conflicts with transcript context')
        out.append(instance)
    write_rows(args.output, out)

def peptide_keys(args):
    seen = {}
    for x in rows(args.instances):
        key = f"{x['interval_id']}|{x['translation']['frame']}"
        seen.setdefault(key, {'peptide_id': key, 'interval_id': x['interval_id'], 'frame': x['translation']['frame'],
                              'sequence': x.get('peptide_sequence', ''), 'instance_ids': []})['instance_ids'].append(x['instance_id'])
    write_rows(args.output, seen.values())

def axis(args):
    name, out = args.axis, []
    for x in rows(args.input):
        item = dict(x)
        # Honest default: a missing reference/tool result is uninformative, never negative.
        item['axes'] = item.get('axes', {})
        item['axes'][name] = typed('uninformative', f'{name} requires a pinned external result; no result supplied')
        out.append(item)
    write_rows(args.output, out)

def merge_axes(args):
    merged = {}
    for path in args.inputs:
        for x in rows(path):
            key = x.get('instance_id', x.get('peptide_id'))
            target = merged.setdefault(key, dict(x))
            target.setdefault('axes', {}).update(x.get('axes', {}))
    write_rows(args.output, merged.values())

def adjudicate(args):
    cluster_by_instance = {}
    if args.clusters:
        for cluster in rows(args.clusters):
            for instance_id in cluster.get('instance_ids', []):
                cluster_by_instance[instance_id] = cluster
    by_interval = defaultdict(list)
    for x in rows(args.instances):
        axes = x.get('axes', {})
        states = {k: v.get('state') for k, v in axes.items()}
        if x.get('admissibility', {}).get('state') == 'unattributable': claim = 'dual-frame-unresolved'
        elif x.get('peptide_uniqueness', {}).get('class') in ('identical', 'substring'): claim = 'fragment-of-known-product'
        elif states.get('regulatory_geometry') == 'positive': claim = 'regulatory-translon'
        elif all(v in ('null', 'uninformative') for v in states.values()): claim = 'characterised-null'
        else: claim = 'candidate-coding-ncORF'
        cluster = cluster_by_instance.get(x.get('instance_id'))
        x['claim'] = {'class': claim, 'licence': 'characterisation only; detection never licenses coding annotation',
                      'axis_states': states,
                      'product_cluster': None if cluster is None else {'sequence_clusters': cluster.get('sequence_clusters', []), 'structural_cluster': cluster.get('structural_cluster')}}
        by_interval[x['interval_id']].append(x)
    result = []
    for iid, items in by_interval.items():
        classes = {x['class'] for x in items}
        result.extend([{**x, 'cross_isoform_coherence': 'coherent' if len(classes) == 1 else 'incoherent'} for x in items])
    write_rows(args.output, result)
    if args.manifest:
        axis_states = defaultdict(lambda: defaultdict(int))
        claims = defaultdict(int)
        for item in result:
            claims[item['claim']['class']] += 1
            for name, axis in item.get('axes', {}).items():
                axis_states[name][axis.get('state', 'missing')] += 1
        with open(args.manifest, 'w') as handle:
            json.dump({'record_type': 'translon-characterisation-manifest', 'instance_count': len(result), 'claim_counts': dict(claims), 'axis_state_counts': {name: dict(values) for name, values in axis_states.items()}, 'coding_wall': 'characterisation only; detection never licenses coding annotation'}, handle, sort_keys=True)

def main():
    p = argparse.ArgumentParser(); sub = p.add_subparsers(required=True)
    s = sub.add_parser('substrate'); s.add_argument('--intervals', required=True); s.add_argument('--gff', required=True); s.add_argument('--genome', required=True); s.add_argument('--instances', required=True); s.add_argument('--unhosted', required=True); s.set_defaults(func=substrate)
    s = sub.add_parser('translation-join'); s.add_argument('--instances', required=True); s.add_argument('--verdicts', required=True); s.add_argument('--output', required=True); s.set_defaults(func=translation_join)
    s = sub.add_parser('peptide-keys'); s.add_argument('--instances', required=True); s.add_argument('--output', required=True); s.set_defaults(func=peptide_keys)
    s = sub.add_parser('axis'); s.add_argument('--axis', required=True); s.add_argument('--input', required=True); s.add_argument('--output', required=True); s.set_defaults(func=axis)
    s = sub.add_parser('merge-axes'); s.add_argument('--inputs', required=True, nargs='+'); s.add_argument('--output', required=True); s.set_defaults(func=merge_axes)
    s = sub.add_parser('adjudicate'); s.add_argument('--instances', required=True); s.add_argument('--clusters'); s.add_argument('--output', required=True); s.add_argument('--manifest'); s.set_defaults(func=adjudicate)
    a = p.parse_args(); a.func(a)
if __name__ == '__main__': main()
