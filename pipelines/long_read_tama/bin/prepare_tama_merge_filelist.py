#!/usr/bin/env python
import argparse
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument("output")
    p.add_argument("beds", nargs="+")
    a = p.parse_args()

    with open(a.output, "w") as handle:
        for bed in sorted(a.beds):
            stem = os.path.splitext(os.path.basename(bed))[0]
            handle.write("{}\t{}\n".format(stem, bed))


if __name__ == "__main__":
    main()