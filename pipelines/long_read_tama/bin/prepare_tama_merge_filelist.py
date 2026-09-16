#!/usr/bin/env python3
import argparse
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("output")
    p.add_argument("beds", nargs="+")
    a = p.parse_args()
    if not a.beds:
        raise SystemExit("No TAMA BED inputs available")
    with open(a.output, "w") as handle:
        for bed in sorted(a.beds):
            handle.write(f"{Path(bed).stem}\t{bed}\n")


if __name__ == "__main__":
    main()
