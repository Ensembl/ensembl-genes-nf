#!/usr/bin/env python
import argparse
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("output")
    p.add_argument("beds", nargs="+")
    a = p.parse_args()

    with open(a.output, "w") as handle:
        for bed in sorted(a.beds):
            handle.write("{}\t{}\n".format(Path(bed).stem, bed))


if __name__ == "__main__":
    main()