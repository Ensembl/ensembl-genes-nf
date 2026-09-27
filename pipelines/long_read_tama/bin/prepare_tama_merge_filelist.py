#!/usr/bin/env python

import argparse
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument("output")
    p.add_argument("beds", nargs="+")
    p.add_argument("--seq-type", default="no_cap")
    p.add_argument("--priority-rank", default="1,1,1")
    args = p.parse_args()

    with open(args.output, "w") as handle:
        for bed in sorted(args.beds):
            filename = os.path.basename(bed)
            stem = os.path.splitext(filename)[0]
            parent = os.path.basename(os.path.dirname(bed))

            # validated_tama.bed is repeated across shards, so use bed01, bed02, ...
            source_id = parent if stem == "validated_tama" and parent else stem

            handle.write("{}\t{}\t{}\t{}\n".format(
                bed, args.seq_type, args.priority_rank, source_id
            ))


if __name__ == "__main__":
    main()
