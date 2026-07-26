Statistics Pipeline
===================

The **Statistics Pipeline** generates comprehensive quality metrics and statistics
for gene annotations and assemblies. This pipeline is essential for validating
annotations, assessing completeness, and providing metadata for Ensembl databases.

Quick Links
-----------

* :doc:`Input Specification <input>`
* :doc:`Output Reference <output>`
* :doc:`Troubleshooting <troubleshooting>`


Workflow Documentation
----------------------

.. list-table::
    :header-rows: 1
    :widths: 25 35 40

   * - Workflow
        - Purpose
        - Key Features
   * - :doc:`BUSCO <workflows/busco>`
        - Assess annotation/assembly completeness
        - Single-copy ortholog presence, protein & genome modes
   * - :doc:`OMArk <workflows/omark>`
        - Proteome quality & contamination screening
        - Consistency checks, lineage validation
   * - :doc:`Ensembl Stats <workflows/ensembl-stats>`
        - Generate database statistics
        - Gene counts, transcript metrics, metakeys

Module Documentation
--------------------

.. toctree::
    :maxdepth: 1

    modules/fetch-genome
    modules/fetch-proteins
    modules/busco-dataset
    modules/busco-genome-lineage
    modules/busco-protein-lineage
    modules/busco-core-metakeys
    modules/omamer-hog
    modules/omark
    modules/run-statistics
    modules/run-ensembl-meta
    modules/populate-db
    modules/db-metadata


Workflow Selection Guide
------------------------

::

Need to assess...
├─ Assembly quality?
│  └─ Use BUSCO (genome mode)
├─ Annotation completeness?
│  └─ Use BUSCO (protein mode)
├─ Contamination?
│  └─ Use OMArk
├─ Database statistics?
│  └─ Use Ensembl Stats
└─ Complete QC?
└─ Use all three workflows
```

Common Use Cases
----------------

Complete Quality Control
^^^^^^^^^^^^^^^^^^^^^^^^

Run all workflows for comprehensive assessment:

.. code-block:: bash

    nextflow run main.nf \
    --csvFile genomes.csv \
    --run_busco_core \
    --busco_mode both \
    --run_omark \
    --run_ensembl_stats \
    --host mysql-server.example.com \
    --user_r ensro \
    --enscode /path/to/ENSCODE \
    --outdir qc_results

.. note::

Provides:
* ✅ Assembly completeness (BUSCO genome)
* ✅ Annotation completeness (BUSCO protein)
* ✅ Contamination screening (OMArk)
* ✅ Database statistics (Ensembl Stats)

Pre-Release Validation
^^^^^^^^^^^^^^^^^^^^^^

Validate before public release:

.. code-block:: bash

    nextflow run main.nf \
    --csvFile release_databases.csv \
    --run_busco_core \
    --busco_mode protein \
    --run_omark \
    --run_ensembl_stats \
    --apply_ensembl_stats \
    --host staging-db.example.com \
    --user ensadmin \
    --password ${DB_PASS} \
    --enscode /path/to/ENSCODE \
    --team genebuild \
    --outdir release_validation


NCBI Assembly Assessment
^^^^^^^^^^^^^^^^^^^^^^^^

Download and assess NCBI assemblies:

.. code-block:: bash

* Create CSV with NCBI assembly accessions
cat > ncbi_assemblies.csv << EOF
dbname,species_id,taxon_id,assembly_accession,assembly_name,taxon_name
gca_001234567_core,1,9606,GCA_001234567.1,ASM123456v1,homo_sapiens
gca_002345678_core,1,10090,GCA_002345678.1,ASM234567v1,mus_musculus
EOF

nextflow run main.nf \
  --csvFile ncbi_assemblies.csv \
  --run_busco_ncbi \
  --outdir ncbi_assessment


Comparative Analysis
^^^^^^^^^^^^^^^^^^^^

Compare quality across multiple species:

.. code-block:: bash

# Vertebrate comparison
cat > vertebrates.csv << EOF
dbname,species_id,taxon_id
homo_sapiens_core_110_38,1,9606
mus_musculus_core_110_39,1,10090
gallus_gallus_core_110_7,1,9031
danio_rerio_core_110_11,1,7955
EOF

nextflow run main.nf \
  --csvFile vertebrates.csv \
  --run_busco_core \
  --busco_mode protein \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro \
  --outdir vertebrate_comparison


Pipeline Architecture
^^^^^^^^^^^^^^^^^^^^
.. code-block:: text

Input CSV
    │
    ├─── BUSCO Analysis
    │    ├─ Protein mode → Annotation completeness
    │    └─ Genome mode → Assembly completeness
    │
    ├─── OMArk Analysis
    │    ├─ Completeness assessment
    │    └─ Contamination detection
    │
    └─── Ensembl Stats
        ├─ Gene/transcript counts
        ├─ Biotype distributions
        └─ Metakey generation
            │
            └─ Apply to database (optional)


Key Features
------------

Flexible Input Options
^^^^^^^^^^^^^^^^^^^^^^

* **Core databases**: Connect to existing Ensembl core databases
* **NCBI assemblies**: Automatically download and assess
* **Mixed sources**: Combine different input types

Multiple Analysis Modes
^^^^^^^^^^^^^^^^^^^^^^^

* **BUSCO**: Protein, genome, or both modes
* **OMArk**: Proteome-based quality with contamination detection
* **Ensembl Stats**: Comprehensive database metrics

Database Integration
^^^^^^^^^^^^^^^^^^^^^^^

*-* **Read-only mode**: Generate statistics without modifying databases
* **Apply mode**: Load statistics and metakeys into databases
* **Validation**: Pre-check before applying changes

Batch Processing
^^^^^^^^^^^^^^^^^^^^^^^

*  Process hundreds of genomes in parallel
*  Automatic resource management
*  Resume capability for interrupted runs

Output Overview
---------------

Directory Structure
^^^^^^^^^^^^^^^^^^^

.. code-block:: text

results/
├── busco/
│   ├── sample1_busco_short_summary.txt
│   ├── sample1_genome_busco_short_summary.txt
│   └── sample1_busco_full_table.tsv
├── omark/
│   └── sample1_omark_proteins_detailed_summary.txt
└── ensembl_stats/
    └── sample1_statistics.json


Result Interpretation
^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :header-rows: 1
    :widths: 30 20 20 20

    * - Metric
        - Good Range
        - Warning Range
        - Action Needed
    * - **BUSCO Complete**
        - >95%
        - 85–95%
        - <85%
    * - **OMArk Consistency**
        - >98%
        - 95–98%
        - <95%
    * - **Gene Count**
        - Expected ±10%
        - Expected ±20%
        - Outside ±20%

Requirements
------------

System Requirements
^^^^^^^^^^^^^^^^^^^

* **Nextflow**: 24.10.3 or higher
* **Java**: 11 or higher
* **Memory**: 32+ GB recommended
* **Storage**: 50+ GB for temporary files

Software Dependencies
^^^^^^^^^^^^^^^^^^^^^

* **BUSCO**: 6.0.0+
* **OMArk**: Latest version
* **Ensembl API**: Release-specific
* **Singularity/Docker**: For containerized workflows

Database Access
^^^^^^^^^^^^^^^
* **MySQL client**: For core database access
* **Read access**: For statistics generation
* **Write access**: For applying metakeys (optional)

Getting Started
---------------

Install Nextflow
^^^^^^^^^^^^^^^^

.. code-block:: bash

curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/


Clone Pipeline
^^^^^^^^^^^^^^

.. code-block:: bash

git clone https://github.com/Ensembl/ensembl-genes.git
cd ensembl-genes/statistics


Prepare Input
^^^^^^^^^^^^^

Create a CSV file with your targets:

.. code-block:: csv
dbname,species_id,taxon_id
homo_sapiens_core_110_38,1,9606


Run Pipeline
^^^^^^^^^^^^

.. code-block:: bash

    nextflow run main.nf \
    --csvFile input.csv \
    --run_busco_core \
    --run_omark \
    --host mysql-server.example.com \
    --user_r ensro \
    --outdir results

Review Results
^^^^^^^^^^^^^^

.. code-block:: bash

# Check BUSCO completeness
cat results/busco/*_short_summary.txt

# Check OMArk consistency
cat results/omark/*_detailed_summary.txt

# Review statistics
cat results/ensembl_stats/*.json


Best Practices
--------------

.. tip::

   **Run All Workflows**
    For production annotations, always run BUSCO, OMArk, and Ensembl Stats together for comprehensive QC.

.. tip::

   **Validate Before Applying**

    Generate and review statistics before using `--apply_*` flags to load data into databases.

.. tip::

   **Use Specific Lineages**
    Choose the most specific BUSCO/OMArk lineage for your organism for best results.

.. tip::

   **Track Over Time**
    Keep statistics outputs in version control to monitor quality trends across releases.

.. tip::

   **Document Exceptions**
    Some species have genuine biological variations (gene losses, duplications) that affect scores—document these.

Common Workflows by Role
------------------------

Annotation Completeness
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

# Complete annotation QC
    nextflow run main.nf \
    --csvFile new_annotations.csv \
    --run_busco_core \
    --busco_mode both \
    --run_omark \
    --run_ensembl_stats \
    --host mysql-server.example.com \
    --user_r ensro \
    --enscode /software/ensembl/ENSCODE


Assembly Completeness
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

# Assembly quality assessment
nextflow run main.nf \
    --csvFile assemblies.csv \
    --run_busco_ncbi \
    --outdir assembly_qc


Ensembl Statistics
^^^^^^^^^^^^^^^^^^

.. code-block:: bash

# Generate and apply statistics
nextflow run main.nf \
    --csvFile release_dbs.csv \
    --run_ensembl_stats \
    --apply_ensembl_stats \
    --host mysql-server.example.com \
    --user ensadmin \
    --password ${DB_PASS} \
    --enscode /software/ensembl/ENSCODE \
    --team genebuild

Module Overview
---------------

The statistics pipeline consists of 13 modules organized into functional categories:

Data Retrieval Modules
^^^^^^^^^^^^^^^^^^^^^^
#. :doc:`fetch-genome <modules/fetch-genome>` – Retrieves genome sequences from Ensembl core databases.
#. :doc:`fetch-proteins <modules/fetch-proteins>` – Extracts protein translations from Ensembl databases.

BUSCO Quality Assessment Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. :doc:`busco-dataset <modules/busco-dataset>` – Downloads appropriate BUSCO lineage datasets.
#. :doc:`busco-genome-lineage <modules/busco-genome-lineage>` – Runs BUSCO assessment on genome sequences.
#. :doc:`busco-protein-lineage <modules/busco-protein-lineage>` – Runs BUSCO assessment on protein translations.
#. :doc:`busco-core-metakeys <modules/busco-core-metakeys>` – Patches BUSCO metadata into core databases.

Orthology Analysis Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^

#. :doc:`omamer-hog <modules/omamer-hog>` – Performs orthology inference using OMAmer.
#. :doc:`omark <modules/omark>` – Quality assessment of protein annotations using OMArk.

Statistics Generation Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. :doc:`run-statistics <modules/run-statistics>` – Generates comprehensive annotation statistics.
#. :doc:`run-ensembl-meta <modules/run-ensembl-meta>` – Generates core database metadata SQL files.

Database Operations Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. :doc:`populate-db <modules/populate-db>` – Executes SQL files to populate databases.
#. :doc:`db-metadata <modules/db-metadata>` – Manages database metadata and versioning.


Pipeline Flow
-------------

.. code-block:: text

1. Data Retrieval
    └─> FETCH_GENOME
    └─> FETCH_PROTEINS

2. Quality Assessment (Parallel)
    ├─> BUSCO_DATASET
    │   ├─> BUSCO_GENOME_LINEAGE
    │   └─> BUSCO_PROTEIN_LINEAGE
    │       └─> BUSCO_CORE_METAKEYS
    │
    └─> OMAMER_HOG
        └─> OMARK

3. Statistics Generation
    ├─> RUN_STATISTICS
    └─> RUN_ENSEMBL_META

4. Database Population
    └─> POPULATE_DB

5. Metadata Management
    └─> DB_METADATA

6. Cleanup
    └─> CLEANING

Module Categories by Function
-----------------------------

Quality Control
^^^^^^^^^^^^^^^

* **BUSCO Genome Lineage** – Assesses genome completeness.
* **BUSCO Protein Lineage** – Assesses proteome completeness.
* **OMArk** – Validates annotation consistency using orthology.

Statistics & Metrics
^^^^^^^^^^^^^^^^^^^^

* **Run Statistics** – Computes gene, transcript and protein counts by biotype.
* **Run Ensembl Meta** – Generates schema and species metadata.

Database Management
^^^^^^^^^^^^^^^^^^^

* **BUSCO Core Metakeys** – Inserts BUSCO metrics into databases.
* **Populate DB** – Executes SQL files for statistics insertion.
* **DB Metadata** – Manages database versioning and tracking.

Resource Management
^^^^^^^^^^^^^^^^^^^

* **Cleaning** – Removes cached intermediate files.

Key Dependencies
----------------

External Tools
^^^^^^^^^^^^^^

* **BUSCO** (v6+) – Genome and proteome completeness assessment.
* **OMAmer** – Orthology inference.
* **OMArk** – Annotation quality assessment.

Ensembl Dependencies
^^^^^^^^^^^^^^^^^^^^

* **Ensembl Perl API** – Database access and manipulation.
* **Ensembl Python libraries** – Metadata generation.
* **Ensembl scripts** – Statistics computation and data extraction.

Databases
^^^^^^^^^

* **Ensembl Core Database** – Primary target for statistics.
* **BUSCO Lineage Datasets** – Reference for completeness assessment.
* **OMAmer HOG Database** – Reference for orthology inference.

Common Parameters
-----------------

Database Connection
^^^^^^^^^^^^^^^^^^^

* ``params.host`` – Database host.
* ``params.port`` – Database port.
* ``params.user`` – Database username.
* ``params.password`` – Database password.

Paths
^^^^^

* ``params.outdir`` – Output directory.
* ``params.cacheDir`` – Cache directory.
* ``params.enscode`` – Path to the Ensembl code repository.

Execution Control
^^^^^^^^^^^^^^^^^

* ``params.files_latency`` – Delay after file operations.
* ``maxForks`` – Maximum number of parallel processes.

Caching Strategy
----------------

Several modules use ``storeDir`` for persistent caching:

* **fetch-genome** – Genome sequences by GCA.
* **fetch-proteins** – Protein translations by GCA.
* **busco-dataset** – BUSCO lineage datasets.
* **omamer-hog** – Orthology assignments by GCA.

This strategy reduces redundant computation and database queries when the
same genomes are processed multiple times.

Conditional Execution
---------------------

Some modules execute conditionally based on parameters:

* **BUSCO_CORE_METAKEYS** – ``params.apply_busco_metakeys``
* **POPULATE_DB** – ``params.apply_ensembl_stats`` or
    ``params.apply_ensembl_beta_metakeys``

Output Structure
----------------

.. code-block:: text

${params.outdir}/
└── ${meta.gca}/
    ├── busco_genome_lineage/
    ├── busco_protein_lineage/
    ├── omark_output/
    ├── core_statistics/
    │   └── *.sql
    └── versions.yml
```

Metadata Requirements
---------------------

* ``gca`` – Genome assembly accession.
* ``dbname`` – Ensembl core database name.
* ``production_name`` – Species production name.
* ``species_id`` – Species identifier.

Documentation Format
--------------------

Each module documentation includes:

* Overview
* Process details
* Inputs
* Outputs
* Parameters
* Script details
* Dependencies
* Notes and best practices

Version Tracking
----------------

All modules generate ``versions.yml`` files containing:

* Tool versions (BUSCO, OMAmer, OMArk)
* Language versions (Python, Perl)
* Database client versions (MySQL)

For More Information
--------------------

* See the individual module documentation.
* Refer to the main pipeline documentation.
* Consult the Ensembl database documentation.
* See the BUSCO, OMAmer and OMArk documentation.

Support
-------

Documentation
^^^^^^^^^^^^^

* Workflow guides
* API reference
* Troubleshooting
* Examples

Getting Help
^^^^^^^^^^^^

* `GitHub Issues <https://github.com/Ensembl/ensembl-genes/issues>`_
* `GitHub Discussions <https://github.com/Ensembl/ensembl-genes/discussions>`_
* Ensembl Genebuild Team

Related Documentation
---------------------

* `Nextflow Documentation <https://www.nextflow.io/docs/latest/>`_
* `BUSCO Documentation <https://busco.ezlab.org/busco_userguide.html>`_
* `OMArk Documentation <https://github.com/DessimozLab/OMArk>`_
* `Ensembl API Documentation <https://www.ensembl.org/info/docs/api/>`_

Citation
--------

If you use this pipeline, please cite:

.. code-block:: text

    Ensembl Genes Statistics Pipeline
    https://github.com/Ensembl/ensembl-genes

Relevant tools:

* **BUSCO** – Manni *et al.* (2021). DOI: ``10.1093/molbev/msab199``
* **OMArk** – Nevers *et al.* (2022). DOI: ``10.1101/2022.11.25.517970``
* **Ensembl** – Cunningham *et al.* (2022). DOI: ``10.1093/nar/gkab1049``
