#!/usr/bin/env python3
"""Stitch exon-level MAF extracts into one strand-correct feature alignment.

Ported in concept from JackCurragh/PhyloCSF-nf, but keyed by the attributed
instance and deliberately limited to an already frame-validated BED12 input.
Missing species are gap-padded in each exon block, preserving codon adjacency.
"""
import argparse
from collections import OrderedDict
from pathlib import Path

COMP = str.maketrans('ACGTNacgtn-', 'TGCANtgcan-')
def rc(s): return s.translate(COMP)[::-1]
def blocks(path):
    block=[]
    for line in Path(path).read_text().splitlines()+['']:
        if line.startswith('s '):
            fields=line.split(); block.append((fields[1].split('.')[0],fields[6]))
        elif not line.strip() and block:
            yield block; block=[]
def main():
    p=argparse.ArgumentParser(); p.add_argument('--mafs',nargs='+',required=True); p.add_argument('--strand',choices=['+','-'],required=True); p.add_argument('--output',required=True); p.add_argument('--species-map'); a=p.parse_args()
    species_map = None
    if a.species_map:
        species_map = dict(line.rstrip().split('\t', 1) for line in Path(a.species_map).read_text().splitlines() if line.strip())
    stitched=OrderedDict()
    for maf in a.mafs:
        for block in blocks(maf):
            width=len(block[0][1]); present={species_map.get(k, k): v for k, v in block if species_map is None or k in species_map}
            if not present:
                continue
            for species in present:
                stitched.setdefault(species, '')
            for species in stitched:
                stitched[species] += present.get(species, '-'*width)
    with open(a.output,'w') as out:
        for species,seq in stitched.items():
            sequence = rc(seq) if a.strand == '-' else seq
            out.write(f'>{species}\n{sequence}\n')
if __name__=='__main__': main()
