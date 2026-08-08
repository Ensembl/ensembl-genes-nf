#!/usr/bin/env python3
"""Materialise and index only the chromosome MAFs required by this run."""
import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fetch(url, destination, expected, algorithm='sha256'):
    last_error = None
    for _ in range(3):
        try:
            with urllib.request.urlopen(url) as source, open(destination, 'wb') as target:
                shutil.copyfileobj(source, target)
            break
        except OSError as error:
            last_error = error
    else:
        raise last_error
    observed = digest(destination, algorithm)
    if observed.lower() != expected.lower():
        destination.unlink(missing_ok=True)
        raise ValueError(f'{algorithm.upper()} mismatch for {url}: expected {expected}, observed {observed}')


def rows(path):
    with open(path) as handle:
        yield from (json.loads(line) for line in handle if line.strip())


def fasta_sizes(path, destination):
    opener = gzip.open if str(path).endswith('.gz') else open
    sizes, name, length = [], None, 0
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    sizes.append((name, length))
                name, length = line[1:].split()[0], 0
            else:
                length += len(line.strip())
    if name is not None:
        sizes.append((name, length))
    if not sizes:
        raise ValueError('reference FASTA contains no sequences')
    destination.write_text(''.join(f'{name}\t{length}\n' for name, length in sizes))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True, help='JSONL: chrom, maf_url and maf_sha256 or maf_md5')
    p.add_argument('--instances', required=True)
    p.add_argument('--genome', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--report', required=True)
    p.add_argument('--maf-index', default='mafIndex')
    a = p.parse_args()
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    required_chroms = {item['chrom'] for item in rows(a.instances)}
    if not required_chroms:
        raise ValueError('cannot prepare MAF reference for an empty instance set')
    chrom_sizes = outdir / 'reference.chrom.sizes'
    fasta_sizes(a.genome, chrom_sizes)
    report, found = [], set()
    for line in Path(a.manifest).read_text().splitlines():
        item = json.loads(line)
        if item.get('chrom') not in required_chroms:
            continue
        required = ('chrom', 'maf_url')
        missing = [key for key in required if not item.get(key)]
        if missing:
            raise ValueError(f'MAF manifest row lacks {missing}')
        checksum, algorithm = (item.get('maf_sha256'), 'sha256') if item.get('maf_sha256') else (item.get('maf_md5'), 'md5')
        if not checksum:
            raise ValueError(f"MAF manifest row for {item['chrom']} has no checksum")
        found.add(item['chrom'])
        maf = outdir / f"{item['chrom']}.maf"
        index = outdir / f"{item['chrom']}.maf.bb"
        archive = outdir / Path(item['maf_url']).name
        fetch(item['maf_url'], archive, checksum, algorithm)
        if archive.suffix == '.gz':
            with gzip.open(archive, 'rb') as source, open(maf, 'wb') as destination:
                shutil.copyfileobj(source, destination)
            archive.unlink()
        else:
            archive.replace(maf)
        subprocess.run([a.maf_index, str(maf), str(index), f'-chromSizes={chrom_sizes}'], check=True)
        report.append({'chrom': item['chrom'], 'maf': maf.name, 'index': index.name,
                       'source_checksum': {'algorithm': algorithm, 'value': checksum},
                       'materialised_sha256': {'maf': digest(maf), 'index': digest(index)}})
    missing_chroms = sorted(required_chroms - found)
    if missing_chroms:
        raise ValueError(f'MAF manifest lacks required chromosomes: {missing_chroms}')
    Path(a.report).write_text(''.join(json.dumps(row, sort_keys=True) + '\n' for row in report))


if __name__ == '__main__':
    main()
