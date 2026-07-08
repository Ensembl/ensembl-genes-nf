#!/usr/bin/env python3

import argparse
import gzip
import heapq
import json
import tempfile
from array import array
from pathlib import Path
from typing import Iterator, Optional, Tuple, List

import numpy as np


UINT32_MAX = np.iinfo(np.uint32).max


def read_sorted_tsv(path: Path) -> Iterator[Tuple[str, int]]:
    """
    Yield coalesced (sequence, count) from a lexically sorted TSV.
    """
    with open(path, "rt") as f:
        previous_seq: Optional[str] = None
        previous_count = 0

        for line_num, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"Malformed line {path}:{line_num}: {line!r}")

            seq = parts[0]
            count = int(parts[1])

            if count < 0:
                raise ValueError(f"Negative count {path}:{line_num}: {count}")

            if previous_seq is not None and seq < previous_seq:
                raise ValueError(
                    f"Input not sorted at {path}:{line_num}: "
                    f"{seq!r} follows {previous_seq!r}"
                )

            if seq == previous_seq:
                previous_count += count
            else:
                if previous_seq is not None:
                    yield previous_seq, previous_count
                previous_seq = seq
                previous_count = count

        if previous_seq is not None:
            yield previous_seq, previous_count


class UInt32Spool:
    """
    Append uint32 values in chunks to a raw binary file.
    Avoids holding enormous Python/NumPy arrays in RAM.
    """

    def __init__(self, path: Path, chunk_size: int = 1_000_000):
        self.path = path
        self.chunk_size = chunk_size
        self.buffer = array("I")
        self.n = 0
        self.fh = open(path, "wb")

    def append(self, value: int):
        if value > UINT32_MAX:
            raise OverflowError(f"uint32 overflow: {value}")
        self.buffer.append(value)
        self.n += 1
        if len(self.buffer) >= self.chunk_size:
            self.flush()

    def flush(self):
        if self.buffer:
            self.buffer.tofile(self.fh)
            self.buffer = array("I")

    def close(self):
        self.flush()
        self.fh.close()

    def memmap(self):
        if self.n == 0:
            return np.asarray([], dtype=np.uint32)
        return np.memmap(self.path, dtype=np.uint32, mode="r", shape=(self.n,))


def kway_merge_counts(tsv_paths: List[Path]) -> Iterator[Tuple[str, List[Tuple[int, int]]]]:
    """
    Merge sorted sample TSVs directly.

    Yields:
        sequence, [(sample_idx, count), ...]
    """
    iterators = []
    heap = []

    for sample_idx, path in enumerate(tsv_paths):
        it = read_sorted_tsv(path)
        iterators.append(it)

        try:
            seq, count = next(it)
            heapq.heappush(heap, (seq, sample_idx, count))
        except StopIteration:
            pass

    while heap:
        seq, sample_idx, count = heapq.heappop(heap)
        counts = [(sample_idx, count)]

        try:
            next_seq, next_count = next(iterators[sample_idx])
            heapq.heappush(heap, (next_seq, sample_idx, next_count))
        except StopIteration:
            pass

        while heap and heap[0][0] == seq:
            _, sample_idx, count = heapq.heappop(heap)
            counts.append((sample_idx, count))

            try:
                next_seq, next_count = next(iterators[sample_idx])
                heapq.heappush(heap, (next_seq, sample_idx, next_count))
            except StopIteration:
                pass

        counts.sort()
        yield seq, counts


def write_csr_npz(path: Path, data, indices, indptr, shape):
    """
    Write a SciPy-compatible CSR .npz without constructing a scipy.sparse object.
    scipy.sparse.load_npz(path) should read this.
    """
    np.savez_compressed(
        path,
        data=data,
        indices=indices,
        indptr=indptr,
        shape=np.asarray(shape, dtype=np.int64),
        format=np.asarray("csr", dtype="|S3"),
    )


def build_matrix(
    study_id: str,
    sample_ids: List[str],
    tsv_paths: List[Path],
    output_dir: Path,
    partition: Optional[str] = None,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = f"{study_id}.{partition}" if partition else study_id
    n_samples = len(sample_ids)

    matrix_path = output_dir / f"{prefix}_matrix.npz"
    sequences_path = output_dir / f"{prefix}_sequences.txt.gz"
    metadata_path = output_dir / f"{prefix}_metadata.json"

    print(f"Building study matrix for {study_id}", flush=True)
    print(f"  Samples: {n_samples}", flush=True)
    print("  Streaming k-way merge", flush=True)

    with tempfile.TemporaryDirectory(prefix=f"{prefix}.csr.") as td:
        td = Path(td)

        data_spool = UInt32Spool(td / "data.u32")
        indices_spool = UInt32Spool(td / "indices.u32")
        indptr_spool = UInt32Spool(td / "indptr.u32")

        n_reads = 0
        nnz = 0
        total_counts = 0

        indptr_spool.append(0)

        with gzip.open(sequences_path, "wt") as seq_out:
            for seq, counts in kway_merge_counts(tsv_paths):
                seq_out.write(seq + "\n")

                for sample_idx, count in counts:
                    if sample_idx >= n_samples:
                        raise ValueError(f"sample_idx out of range: {sample_idx}")

                    indices_spool.append(sample_idx)
                    data_spool.append(count)
                    nnz += 1
                    total_counts += count

                n_reads += 1
                indptr_spool.append(nnz)

                if n_reads % 1_000_000 == 0:
                    print(f"  Processed {n_reads:,} unique reads; nnz={nnz:,}", flush=True)

        data_spool.close()
        indices_spool.close()
        indptr_spool.close()

        print(f"  Final rows: {n_reads:,}", flush=True)
        print(f"  Final nnz: {nnz:,}", flush=True)
        print("  Writing CSR npz", flush=True)

        write_csr_npz(
            matrix_path,
            data=data_spool.memmap(),
            indices=indices_spool.memmap(),
            indptr=indptr_spool.memmap(),
            shape=(n_reads, n_samples),
        )

    metadata = {
        "study_id": study_id,
        "partition": partition,
        "samples": sample_ids,
        "n_reads": n_reads,
        "n_samples": n_samples,
        "nnz": nnz,
        "total_counts": int(total_counts),
        "matrix_density": nnz / (n_reads * n_samples) if n_reads and n_samples else 0,
        "row_identity": "row_id is zero-based line number in sequences.txt.gz",
        "matrix_format": "scipy CSR npz",
        "vocab": "not written; use sequences file as row-id mapping",
    }

    with open(metadata_path, "wt") as f:
        json.dump(metadata, f, indent=2)

    print("Outputs:", flush=True)
    print(f"  Matrix:    {matrix_path}", flush=True)
    print(f"  Sequences: {sequences_path}", flush=True)
    print(f"  Metadata:  {metadata_path}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv_files", nargs="+", type=Path)
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument("-s", "--study-id", required=True)
    parser.add_argument("--sample-ids")
    parser.add_argument("--partition")

    args = parser.parse_args()

    if args.sample_ids:
        sample_ids = args.sample_ids.split(",")
        if len(sample_ids) != len(args.tsv_files):
            raise SystemExit(
                f"{len(sample_ids)} sample IDs for {len(args.tsv_files)} TSV files"
            )
    else:
        sample_ids = []
        for path in args.tsv_files:
            stem = path.name
            if stem.endswith(".tsv"):
                stem = stem[:-4]
            if args.partition and stem.endswith(f".{args.partition}"):
                stem = stem[: -len(f".{args.partition}")]
            sample_ids.append(stem)

    build_matrix(
        study_id=args.study_id,
        sample_ids=sample_ids,
        tsv_paths=args.tsv_files,
        output_dir=args.output_dir,
        partition=args.partition,
    )


if __name__ == "__main__":
    main()
