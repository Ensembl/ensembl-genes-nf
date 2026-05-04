#!/usr/bin/env python3
import argparse
import contextlib
import gzip
import heapq
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple


# All possible dinucleotide prefixes (sorted for consistent ordering)
DINUCLEOTIDES = sorted([f"{a}{b}" for a in "ACGT" for b in "ACGT"])
PARTITIONS = DINUCLEOTIDES + ["NN"]
COUNT_PATTERN = re.compile(r"_x(\d+)(?:\s|$)")


def get_dinucleotide(seq: str) -> str:
    """Get dinucleotide prefix, defaulting to 'NN' for short/invalid sequences."""
    if len(seq) >= 2:
        prefix = seq[:2].upper()
        if prefix in DINUCLEOTIDES:
            return prefix
    return "NN"


def open_text(path: Path):
    """Open plain or gzipped text input."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_count(header: str) -> int:
    """Parse collapsed FASTA count from a header such as >read_x1234."""
    match = COUNT_PATTERN.search(header)
    if match:
        return int(match.group(1))

    print(f"Warning: Could not parse count from header: {header}", file=sys.stderr)
    return 1


def iter_collapsed_fasta(fasta_path: Path) -> Iterator[Tuple[str, int]]:
    """
    Stream collapsed FASTA records as (sequence, count).

    Multi-line FASTA records are accepted. Duplicate sequences are deliberately
    handled later by the sorter, so this iterator never retains the whole sample.
    """
    count: Optional[int] = None
    sequence_parts: List[str] = []

    def flush_record() -> Optional[Tuple[str, int]]:
        if count is None or not sequence_parts:
            return None
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            return None
        return sequence, count

    with open_text(fasta_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                record = flush_record()
                if record is not None:
                    yield record
                count = parse_count(line)
                sequence_parts = []
            else:
                sequence_parts.append(line)

    record = flush_record()
    if record is not None:
        yield record


def read_sorted_tsv(path: Path) -> Iterator[Tuple[str, int]]:
    """Read a sorted temporary TSV chunk."""
    with open(path, "r") as handle:
        for line_num, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            try:
                sequence, count = line.split("\t", 1)
                yield sequence, int(count)
            except ValueError as exc:
                raise ValueError(f"Malformed TSV line in {path}:{line_num}: {line}") from exc


def write_atomic_sorted_merge(chunk_paths: List[Path], output_path: Path) -> int:
    """
    Merge sorted chunk files into output_path and coalesce duplicate sequences.

    Returns the number of unique output sequences.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_handle = tempfile.NamedTemporaryFile(
        "w",
        delete=False,
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
    )
    tmp_output = Path(output_handle.name)
    unique_count = 0

    try:
        with output_handle, contextlib.ExitStack() as stack:
            iterators = [stack.enter_context(open(path, "r")) for path in chunk_paths]
            parsed_iterators = (
                (tuple(line.rstrip("\n").split("\t", 1)) for line in handle if line.strip())
                for handle in iterators
            )

            current_seq: Optional[str] = None
            current_count = 0

            for sequence, count_text in heapq.merge(*parsed_iterators, key=lambda item: item[0]):
                count = int(count_text)
                if current_seq is None:
                    current_seq = sequence
                    current_count = count
                elif sequence == current_seq:
                    current_count += count
                else:
                    output_handle.write(f"{current_seq}\t{current_count}\n")
                    unique_count += 1
                    current_seq = sequence
                    current_count = count

            if current_seq is not None:
                output_handle.write(f"{current_seq}\t{current_count}\n")
                unique_count += 1

        os.replace(tmp_output, output_path)
        return unique_count
    except Exception:
        try:
            tmp_output.unlink()
        except FileNotFoundError:
            pass
        raise


class SpillSorter:
    """Memory-bounded sequence/count sorter backed by temporary sorted chunks."""

    def __init__(
        self,
        output_path: Path,
        temp_dir: Path,
        max_records_in_memory: int,
        merge_fan_in: int,
    ):
        self.output_path = output_path
        self.temp_dir = temp_dir
        self.max_records_in_memory = max(1, max_records_in_memory)
        self.merge_fan_in = max(2, merge_fan_in)
        self.buffer: Dict[str, int] = {}
        self.chunk_paths: List[Path] = []
        self.total_records = 0
        self.total_counts = 0
        self.unique_records = 0

    def add(self, sequence: str, count: int):
        self.total_records += 1
        self.total_counts += count
        self.buffer[sequence] = self.buffer.get(sequence, 0) + count
        if len(self.buffer) >= self.max_records_in_memory:
            self.flush()

    def flush(self):
        if not self.buffer:
            return

        chunk = tempfile.NamedTemporaryFile(
            "w",
            delete=False,
            dir=self.temp_dir,
            prefix="collapsed_to_tsv.",
            suffix=".chunk.tsv",
        )
        chunk_path = Path(chunk.name)
        with chunk:
            for sequence, count in sorted(self.buffer.items()):
                chunk.write(f"{sequence}\t{count}\n")

        self.chunk_paths.append(chunk_path)
        self.buffer.clear()

    def _compact_chunks(self) -> List[Path]:
        chunks = self.chunk_paths
        while len(chunks) > self.merge_fan_in:
            next_round: List[Path] = []
            for start in range(0, len(chunks), self.merge_fan_in):
                group = chunks[start:start + self.merge_fan_in]
                merged = tempfile.NamedTemporaryFile(
                    "w",
                    delete=False,
                    dir=self.temp_dir,
                    prefix="collapsed_to_tsv.merge.",
                    suffix=".chunk.tsv",
                )
                merged_path = Path(merged.name)
                merged.close()
                write_atomic_sorted_merge(group, merged_path)
                next_round.append(merged_path)
                for chunk in group:
                    try:
                        chunk.unlink()
                    except FileNotFoundError:
                        pass
            chunks = next_round
        return chunks

    def finish(self) -> int:
        self.flush()

        if not self.chunk_paths:
            self.output_path.parent.mkdir(parents=True, exist_ok=True)
            self.output_path.write_text("")
            self.unique_records = 0
            return 0

        chunks = self._compact_chunks()
        try:
            self.unique_records = write_atomic_sorted_merge(chunks, self.output_path)
            return self.unique_records
        finally:
            for chunk in chunks:
                try:
                    chunk.unlink()
                except FileNotFoundError:
                    pass


def main():
    parser = argparse.ArgumentParser(
        description='Convert collapsed FASTA to sorted TSV(s)'
    )
    parser.add_argument(
        'input',
        type=Path,
        help='Input collapsed FASTA file (.fa or .fa.gz)'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        help='Output TSV file (default: {input_stem}.tsv). Ignored with --partition.'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('.'),
        help='Output directory for partition files (default: current directory)'
    )
    parser.add_argument(
        '--sample-id',
        type=str,
        help='Sample ID to use in output filename (default: derived from input)'
    )
    parser.add_argument(
        '--partition',
        action='store_true',
        help='Output partition files by dinucleotide prefix (AA, AC, ..., TT, NN)'
    )
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=500_000,
        help='Maximum unique sequences to sort in memory per chunk (default: 500000)'
    )
    parser.add_argument(
        '--merge-fan-in',
        type=int,
        default=256,
        help='Maximum temporary chunk files to merge at once (default: 256)'
    )
    parser.add_argument(
        '--temp-dir',
        type=Path,
        default=None,
        help='Directory for temporary sort chunks (default: output directory)'
    )

    args = parser.parse_args()

    # Determine sample ID
    if args.sample_id:
        sample_id = args.sample_id
    else:
        stem = args.input.name
        if stem.endswith('.gz'):
            stem = stem[:-3]
        if stem.endswith('.fa') or stem.endswith('.fasta'):
            stem = stem.rsplit('.', 1)[0]
        sample_id = stem

    print(f"Reading: {args.input}", file=sys.stderr)
    temp_dir = args.temp_dir or args.output_dir
    temp_dir.mkdir(parents=True, exist_ok=True)

    if args.partition:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Writing {len(PARTITIONS)} partition files to: {args.output_dir}", file=sys.stderr)
        per_partition_chunk_size = max(1, args.chunk_size // len(PARTITIONS))
        sorters = {
            partition: SpillSorter(
                args.output_dir / f"{sample_id}.{partition}.tsv",
                temp_dir,
                per_partition_chunk_size,
                args.merge_fan_in,
            )
            for partition in PARTITIONS
        }

        for sequence, count in iter_collapsed_fasta(args.input):
            sorters[get_dinucleotide(sequence)].add(sequence, count)

        partition_sizes = []
        total_records = 0
        total_counts = 0
        total_unique = 0
        for partition in PARTITIONS:
            sorter = sorters[partition]
            unique_count = sorter.finish()
            partition_sizes.append(unique_count)
            total_records += sorter.total_records
            total_counts += sorter.total_counts
            total_unique += unique_count

        print(f"Read {total_records:,} FASTA records", file=sys.stderr)
        print(f"Found {total_unique:,} unique sequences", file=sys.stderr)
        print(f"Partition sizes: min={min(partition_sizes):,}, max={max(partition_sizes):,}, "
              f"mean={sum(partition_sizes)/len(partition_sizes):,.0f}", file=sys.stderr)
    else:
        if args.output:
            output_path = args.output
        else:
            output_path = args.output_dir / f"{sample_id}.tsv"

        print(f"Writing: {output_path}", file=sys.stderr)
        sorter = SpillSorter(output_path, temp_dir, args.chunk_size, args.merge_fan_in)
        for sequence, count in iter_collapsed_fasta(args.input):
            sorter.add(sequence, count)
        unique_count = sorter.finish()

        print(f"Read {sorter.total_records:,} FASTA records", file=sys.stderr)
        print(f"Found {unique_count:,} unique sequences", file=sys.stderr)
        total_counts = sorter.total_counts
        total_unique = unique_count

    # Summary stats
    print(f"Total read counts: {total_counts:,}", file=sys.stderr)
    if total_unique > 0:
        print(f"Compression ratio: {total_counts / total_unique:.1f}x", file=sys.stderr)
    else:
        print("Compression ratio: N/A (no sequences)", file=sys.stderr)


if __name__ == '__main__':
    main()
