#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', nargs='+', required=True)
    parser.add_argument('--orthology', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    atlas_assemblies = {
        path.name.removeprefix("atlas_")
        for path in (Path(item) for item in args.atlas)
        if path.is_dir()
    }
    with open(args.output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['species', 'assembly', 'ortholog_group', 'status'])
        with open(args.orthology) as orthology:
            for row in csv.DictReader(orthology, delimiter='\t'):
                assembly = row.get('assembly', '')
                status = 'available' if assembly in atlas_assemblies else 'missing_atlas'
                writer.writerow([row.get('species', ''), assembly, row.get('ortholog_group', ''), status])


if __name__ == '__main__':
    main()
