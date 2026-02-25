#!/bin/bash
set -euo pipefail

# Test script for unique reads merger integration
# This script tests both per-run and progressive modes

echo "=== Testing Unique Reads Merger ==="

# Create test directory
TEST_DIR="test_unique_reads_output"
mkdir -p "$TEST_DIR"

# Find some test FASTA files
echo "Finding test files..."
TEST_FILES=$(find work -name "*.collapsed.fa.gz" 2>/dev/null | head -3)

if [ -z "$TEST_FILES" ]; then
    echo "No test files found. Please run the pipeline first to generate some collapsed FASTA files."
    exit 1
fi

# Count how many test files we have
NUM_FILES=$(echo "$TEST_FILES" | wc -l | tr -d ' ')
echo "Found $NUM_FILES test files"

# Create a file list
echo "$TEST_FILES" > "$TEST_DIR/file_list.txt"

echo ""
echo "=== Test 1: Per-Run Mode ==="
bin/sorted_read_index_merger.py \
    --file-list "$TEST_DIR/file_list.txt" \
    --output-dir "$TEST_DIR/per_run" \
    --batch-size 2 \
    --n-workers 2 \
    --output-suffix "_test_run"

# Check outputs
if [ -f "$TEST_DIR/per_run/unique_reads_test_run.fa" ] && \
   [ -f "$TEST_DIR/per_run/count_matrix_test_run.parquet" ] && \
   [ -f "$TEST_DIR/per_run/read_mapping_test_run.json" ]; then
    echo "✓ Per-run mode: All output files created successfully"
else
    echo "✗ Per-run mode: Missing output files"
    exit 1
fi

echo ""
echo "=== Test 2: Progressive Mode (Initial) ==="
bin/sorted_read_index_merger.py \
    --file-list "$TEST_DIR/file_list.txt" \
    --output-dir "$TEST_DIR/progressive_1" \
    --batch-size 2 \
    --n-workers 2

if [ -f "$TEST_DIR/progressive_1/unique_reads.fa" ] && \
   [ -f "$TEST_DIR/progressive_1/count_matrix.parquet" ] && \
   [ -f "$TEST_DIR/progressive_1/read_mapping.json" ]; then
    echo "✓ Progressive mode (initial): All output files created successfully"
else
    echo "✗ Progressive mode (initial): Missing output files"
    exit 1
fi

echo ""
echo "=== Test 3: Progressive Mode (Update) ==="
# For this test, we'll use the same files again to simulate adding more samples
bin/sorted_read_index_merger.py \
    --file-list "$TEST_DIR/file_list.txt" \
    --output-dir "$TEST_DIR/progressive_2" \
    --batch-size 2 \
    --n-workers 2 \
    --previous-index "$TEST_DIR/progressive_1"

if [ -f "$TEST_DIR/progressive_2/unique_reads.fa" ] && \
   [ -f "$TEST_DIR/progressive_2/count_matrix.parquet" ] && \
   [ -f "$TEST_DIR/progressive_2/read_mapping.json" ]; then
    echo "✓ Progressive mode (update): All output files created successfully"
else
    echo "✗ Progressive mode (update): Missing output files"
    exit 1
fi

echo ""
echo "=== Verification ==="

# Check that progressive mode has more samples in the count matrix
echo "Checking count matrix dimensions..."

python3 << 'EOF'
import polars as pl
import json

# Load matrices
df1 = pl.read_parquet("test_unique_reads_output/progressive_1/count_matrix.parquet")
df2 = pl.read_parquet("test_unique_reads_output/progressive_2/count_matrix.parquet")

print(f"Progressive run 1: {len(df1)} unique reads, {len(df1.columns)-1} samples")
print(f"Progressive run 2: {len(df2)} unique reads, {len(df2.columns)-1} samples")

# In this test we used the same files twice, so:
# - Unique reads should be the same (same sequences)
# - Samples should be double (same samples processed twice)
assert len(df2.columns) == len(df1.columns) + len(df1.columns) - 1, \
    f"Expected {len(df1.columns) + len(df1.columns) - 1} columns, got {len(df2.columns)}"

print("✓ Progressive mode correctly merged with previous index")

# Load and check summary
with open("test_unique_reads_output/per_run/processing_summary_test_run.json") as f:
    summary = json.load(f)
    print(f"\nProcessing summary:")
    print(f"  Files processed: {summary['processing_summary']['total_files_processed']}")
    print(f"  Total reads: {summary['processing_summary']['total_reads_processed']:,}")
    print(f"  Unique reads: {summary['processing_summary']['total_unique_reads']:,}")
    print(f"  Compression: {summary['processing_summary']['compression_ratio']:.1f}x")

print("\n✓ All tests passed!")
EOF

echo ""
echo "=== Test Summary ==="
echo "✓ Per-run mode working"
echo "✓ Progressive mode (initial) working"
echo "✓ Progressive mode (update) working"
echo "✓ Output files verified"
echo ""
echo "Test outputs saved in: $TEST_DIR/"
