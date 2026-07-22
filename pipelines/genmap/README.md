# GenMap mappability pipeline

This pipeline calculates reference-genome mappability with [GenMap](https://github.com/cpockrandt/genmap). It uses the pinned `nf-core/genmap/index` and `nf-core/genmap/map` modules and produces text, WIG, bedGraph, and CSV outputs.

The workflow derives a BED interval covering every reference sequence so the GenMap module computes mappability across the complete genome.

The defaults are 28-mers and zero mismatches. GenMap 1.3.0 is supplied by the modules through Bioconda or a container.

## Inputs

Provide one of:

* `--species`: Ensembl species name, for example `homo_sapiens`. The URL is resolved to the species' Ensembl `dna.toplevel.fa.gz` file.
* `--reference_url`: a direct FASTA or FASTA.GZ URL. This takes precedence over Ensembl resolution.
* `--reference`: an existing local FASTA or FASTA.GZ file.

`--release` defaults to `current`; set it to a numbered Ensembl release such as `114` when reproducibility requires a fixed release.

For very large genomes, GenMap's memory-saving index construction can be enabled with `--index_args '-S 20'`. This reduces memory use at the cost of slower indexing.

Human-genome mapping requests 128 GB in the pipeline process configuration. For a site-specific override, provide an additional Nextflow config rather than a pipeline parameter:

```groovy
process {
    withName: 'GENMAP_MAP' {
        memory = 256.GB
    }
}
```

Run with `-c /path/to/site-resources.config`.

## Usage

From this directory, using a container profile:

```bash
nextflow run main.nf -profile docker \
  --species homo_sapiens \
  --release 114 \
  --outdir results/homo_sapiens
```

With a direct URL:

```bash
nextflow run main.nf -profile docker \
  --reference_url 'https://example.org/reference.fa.gz' \
  --outdir results/reference
```

To change the mappability definition:

```bash
nextflow run main.nf -profile docker \
  --species homo_sapiens \
  --kmer 28 \
  --mismatches 0 \
  --outdir results/homo_sapiens
```

Outputs are written beneath `--outdir`:

* `reference/reference.fa` and `reference/reference_source.txt`
* `index/` containing the GenMap index
* `mappability/` containing `.txt`, `.wig`, `.bedgraph`, and `.csv` files

Use `-resume` to reuse a completed download and index.
