# HG002 long-RNA comparison dataset

This is a research inventory derived from Huang, Li and HPRC, “SNP calling,
haplotype phasing and allele-specific analysis with long RNA-seq reads,”
Nature Methods (published 30 March 2026), and its supplementary dataset table.
The paper evaluates PacBio MAS-Seq/Kinnex, PacBio Iso-Seq, ONT cDNA and ONT
direct-RNA data across HG002, HG004 and HG005.

`sample_sheet.tsv` is intentionally an inventory sheet, not an approved
production manifest. The long-read pipeline requires accession-shaped IDs,
selected artifact checksums, exact artifact basenames and a classification
report hash before it will accept `--approved_manifest`. Several paper rows
are exposed through directory listings or BAM URLs rather than stable ENA run
accessions, so those fields must be resolved before production use.

The HG002a MAS/Kinnex row is handled as an explicit reviewed direct-source
exception. Its PacBio directory exposes `segmented.bam` without a published
MD5 or archive accession, so the approved manifest records
`source_type=PACBIO_DIRECT_URI`, a review reason, and
`checksum_policy=record_on_acquisition`. The acquisition task records the
computed MD5 in its checksum report after download; this exception is visible
and is not a validator bypass.

Recommended first comparison:

1. `HG002a_MAS_cDNA` — the PacBio Kinnex/MAS-Seq dataset.
2. `HG002a_ISO_cDNA` — PacBio Iso-Seq counterpart, once its exact GIAB path is
   resolved.
3. `HG002a_ONT_cDNA` — ONT cDNA, after converting the cited BAM to canonical
   FASTQ or adding an approved BAM acquisition route.

Keep direct-RNA rows separate from cDNA rows: they are a useful technology
comparison, but not an interchangeable input because their molecule type and
strand/base-modification properties differ.

The PacBio Kinnex source directory is:

`https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-Revio-HG002-1/`

The paper's complete raw-read and sample descriptions are available from the
[Nature Methods article](https://www.nature.com/articles/s41592-026-03045-6)
and its [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13041723/).
