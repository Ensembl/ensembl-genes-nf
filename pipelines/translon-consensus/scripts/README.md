# Tool-Specific BED12 Conversion Scripts

These scripts convert tool-specific output formats to proper BED12 format for consensus analysis.

## Problem

Different ORF prediction tools output in different formats:
- **RiboTIE**: GTF format with CDS features per exon
- **PRICE**: BED12 (but may have formatting issues)
- **iRibo**: GFF3 format with CDS features per exon
- **ORFQuant**: GFF3 format with CDS features per exon

The consensus pipeline requires proper BED12 format with one entry per ORF (not per exon).

## Scripts

### Auto-Detection (Recommended)

```bash
convert_to_bed12.py input.bed output.bed12
```

Automatically detects the tool format and converts to BED12.

### Tool-Specific Converters

If auto-detection fails, use the tool-specific converter:

#### RiboTIE (GTF → BED12)
```bash
ribotie_to_bed12.py input.gtf output.bed12
```

**Input format:**
```
chr7    RiboTIE CDS     92134173    92134346    .    -    0    gene_id "ENSG..."; ORF_id "ENST..."; ...
chr7    RiboTIE CDS     92131774    92131872    .    -    0    gene_id "ENSG..."; ORF_id "ENST..."; ...
```

**Features:**
- Groups CDS features by `ORF_id`
- Extracts gene_name, transcript_id, ORF_type
- Converts GTF coordinates (1-based) to BED (0-based)
- Creates single BED12 entry per ORF with block structure

#### PRICE (BED12 → BED12)
```bash
price_to_bed12.py input.bed output.bed12
```

**Input format:**
```
1    145842405    145842420    ENST00000334163.4_Trunc_0    80    +    0    0    0    1    15,    0,
```

**Fixes:**
- Validates block count matches actual blocks
- Removes trailing commas (if missing)
- Ensures coordinates are valid
- Fixes thickStart/thickEnd bounds

#### iRibo (GFF3 → BED12)
```bash
iribo_to_bed12.py input.gff3 output.bed12
```

**Input format:**
```
chr1    iRibo    CDS    16529    16765    .    -    .    ID=candidate_orf964261
chr1    iRibo    CDS    16853    17055    .    -    .    ID=candidate_orf964261
```

**Features:**
- Groups CDS features by `ID` attribute
- Converts GFF3 coordinates (1-based) to BED (0-based)
- Creates single BED12 entry per ORF

#### ORFQuant (GFF3 → BED12)
```bash
orfquant_to_bed12.py input.gff3 output.bed12
```

**Input format:**
```
chr1    ORFQuant    CDS    841492    841501    0    +    0    ENST00000670780.1_568_576
chr1    ORFQuant    CDS    1020172   1020373   0    +    0    ENST00000620552.4_51_6254
```

**Features:**
- Groups CDS features by 9th column (ORF identifier)
- Converts GFF3 coordinates (1-based) to BED (0-based)
- Creates single BED12 entry per ORF

## BED12 Format

Output format is standard BED12:

```
chr    start    end    name    score    strand    thickStart    thickEnd    itemRgb    blockCount    blockSizes    blockStarts
```

**Fields:**
- `chr`, `start`, `end`: Genomic location (0-based)
- `name`: ORF identifier (tool-specific)
- `score`: 0-1000 (scaled from tool scores if available)
- `strand`: +/-
- `thickStart`, `thickEnd`: CDS boundaries (same as start/end for ORFs)
- `itemRgb`: RGB color (default: 0,0,0)
- `blockCount`: Number of exons
- `blockSizes`: Comma-separated exon sizes
- `blockStarts`: Comma-separated exon starts (relative to `start`)

## Batch Conversion

Convert all files in a directory:

```bash
for file in input_dir/*.bed; do
    output=$(basename "$file" .bed).bed12
    convert_to_bed12.py "$file" "output_dir/$output"
done
```

## Integration with Pipeline

These scripts should be run **before** the consensus pipeline:

```bash
# 1. Convert all tool outputs to BED12
for tool_dir in RiboTIE PRICE iRibo ORFQuant; do
    for bed in results/$tool_dir/*.bed; do
        base=$(basename "$bed" .bed)
        convert_to_bed12.py "$bed" "results_bed12/$tool_dir/${base}.bed12"
    done
done

# 2. Run consensus pipeline on converted files
nextflow run pipelines/translon-consensus \
    --bed_results_dir results_bed12 \
    --gencode_gtf gencode.v45.annotation.gtf.gz \
    --outdir consensus_output
```

## Testing

Test on example data:

```bash
# Test each converter
cd pipelines/translon-consensus/scripts

# Test auto-detection
./convert_to_bed12.py ../../../test_data/ribotie_sample.bed test_output.bed12

# Verify output
head test_output.bed12
```

## Troubleshooting

**Error: "Could not detect file format"**
- Use `--tool` flag to force a specific converter
- Check that input file has correct structure

**Error: "Block count mismatch"**
- Common with PRICE outputs
- Script will attempt to fix automatically
- Check warnings in stderr

**Error: "Invalid coordinates"**
- Entries with start >= end will be skipped
- Check input file for corruption

## Notes

- All converters handle **multi-exon ORFs** (BED12 blocks)
- Coordinate conversion is handled automatically (1-based → 0-based)
- Duplicate entries with same genomic location will be merged by the consensus pipeline
- Feature names are preserved from input where possible
