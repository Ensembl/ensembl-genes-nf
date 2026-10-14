# BUSCO_LINEAGE Module

## Overview

`BUSCO_LINEAGE` runs BUSCO completeness assessment on either genome assemblies or protein annotations using lineage-specific BUSCO datasets.

The analysis mode is selected dynamically through:

```groovy
meta.busco_mode
```

Supported modes:

| Mode      | Input         | Purpose                        |
| --------- | ------------- | ------------------------------ |
| `genome`  | Genome FASTA  | Assess assembly completeness   |
| `protein` | Protein FASTA | Assess annotation completeness |

**Module Location**

```text
pipelines/statistics/modules/busco_lineage.nf
```

---

## Functionality

The module:

1. Runs BUSCO using a lineage-specific dataset.
2. Measures complete, duplicated, fragmented and missing BUSCO genes.
3. Produces standardized completeness metrics.
4. Operates in offline mode using pre-downloaded datasets.
5. Publishes summary reports and version information.

---

## Inputs

### Channel Input

```groovy
tuple val(meta), path(fasta_file)
```

### Required Metadata

```groovy
[
    gca: String,
    busco_dataset: String,
    busco_mode: String   // genome | protein
]
```

### Input File

| Mode    | Expected Input |
| ------- | -------------- |
| genome  | Genomic FASTA  |
| protein | Protein FASTA  |

---

## Outputs

### Published Files

Output directory:

```text
${params.outdir}/${meta.gca}/
```

Generated files:

```text
busco_<mode>/
├── short_summary.*
├── batch_summary.txt
└── logs/
```

Version file:

```text
versions_busco_<mode>.yml
```

### Output Channels

| Channel                | Description                 |
| ---------------------- | --------------------------- |
| `busco_lineage_output` | BUSCO summary files         |
| `busco_full_output`    | Full BUSCO output directory |
| `versions_file`        | Software versions           |

---

## BUSCO Execution

The module automatically selects the BUSCO mode:

```groovy
def busco_mode_arg =
    meta.busco_mode == 'protein'
        ? 'proteins'
        : 'genome'
```

Core command:

```bash
busco \
    -i ${fasta_file} \
    --mode ${busco_mode_arg} \
    -l ${meta.busco_dataset} \
    -c ${task.cpus} \
    --offline \
    --download_path ${params.download_path}
```

---

## Genome vs Protein Mode

| Feature         | Genome            | Protein          |
| --------------- | ----------------- | ---------------- |
| Input           | DNA FASTA         | Protein FASTA    |
| Gene prediction | Yes               | No               |
| Runtime         | Slower            | Faster           |
| Purpose         | Assembly QC       | Annotation QC    |
| Typical use     | Before annotation | After annotation |

### Genome Mode

BUSCO predicts genes from the assembly before searching for orthologs.

Tools used may include:

* Metaeuk
* Augustus
* Prodigal
* HMMER

Measures assembly completeness.

### Protein Mode

BUSCO searches translated proteins directly.

Tools used:

* HMMER

Measures annotation completeness.

---

## BUSCO Metrics

BUSCO reports:

| Metric | Description                 |
| ------ | --------------------------- |
| C      | Complete BUSCOs             |
| S      | Complete single-copy BUSCOs |
| D      | Complete duplicated BUSCOs  |
| F      | Fragmented BUSCOs           |
| M      | Missing BUSCOs              |

Example:

```text
C:97.2%[S:95.8%,D:1.4%],F:1.3%,M:1.5%
```

---

## Interpretation

### High Quality

```text
C > 95%
M < 2%
```

Indicates a highly complete assembly or annotation.

### Moderate Quality

```text
C = 85–95%
```

Generally suitable for downstream analyses.

### Poor Quality

```text
C < 85%
```

Suggests incomplete assembly or annotation issues.

---

## Comparing Genome and Protein BUSCO

### Good Assembly, Poor Annotation

```text
Genome : 98%
Protein: 82%
```

Interpretation:

* Assembly is complete.
* Gene annotation missed many expected genes.

### Poor Assembly, Better Annotation

```text
Genome : 75%
Protein: 88%
```

Interpretation:

* Assembly quality is limiting.
* Annotation metrics may be misleading.

For robust assessment, both modes should be evaluated.

---

## Process Configuration

```groovy
label 'busco'
maxForks 10
afterScript "sleep ${params.files_latency}"
```

Typical resources:

| CPUs | Memory | Runtime       |
| ---- | ------ | ------------- |
| 8    | 32 GB  | 30 min – 24 h |

---

## Required Parameters

| Parameter              | Description                 |
| ---------------------- | --------------------------- |
| `params.download_path` | BUSCO datasets location     |
| `params.outdir`        | Published results directory |
| `params.files_latency` | File system sync delay      |

---

## Example Workflow

```groovy
FETCH_GENOME
    -> BUSCO_LINEAGE (genome)

FETCH_PROTEINS
    -> BUSCO_LINEAGE (protein)
```

Example metadata:

```groovy
[
    gca: 'GCA_000001405.29',
    busco_dataset: 'primates_odb12',
    busco_mode: 'genome'
]
```

or

```groovy
[
    gca: 'GCA_000001405.29',
    busco_dataset: 'primates_odb12',
    busco_mode: 'protein'
]
```

---

## Version Tracking

Example:

```yaml
"BUSCO_LINEAGE:genome":
  busco: v6.0.0_cv1

"BUSCO_LINEAGE:protein":
  busco: v6.0.0_cv1
```

---

## Related Modules

* [busco-dataset](busco-dataset.md)
* [fetch-genome](fetch-genome.md)
* [fetch-proteins](fetch-proteins.md)

---

**Maintained By:** Ensembl Genes Team
