-- get recent annotations
SELECT
  DISTINCT
  -- genome_annotation.type AS 'genome_annotation.type',
  -- genome_annotation.value AS 'genome_annotation.value',
  -- genome_annotation.genome_database_id AS 'genome_annotation.genome_database_id',
  -- genome_annotation.genome_id AS 'genome_annotation.genome_id',
  -- genome.data_release_id AS 'genome.data_release_id',
  -- data_release.ensembl_genomes_version AS 'data_release.ensembl_genomes_version',
  genome_database.dbname AS 'genome_database.dbname'
FROM genome_annotation
INNER JOIN genome_database
  ON genome_annotation.genome_database_id = genome_database.genome_database_id
INNER JOIN genome
  ON genome_annotation.genome_id = genome.genome_id
INNER JOIN data_release_database
  ON genome.data_release_id = data_release_database.data_release_id
INNER JOIN data_release
  ON data_release_database.data_release_id = data_release.data_release_id
WHERE genome_annotation.type = 'genebuild_method'
  -- select the two most recent released Rapid Release versions
  AND data_release.ensembl_genomes_version BETWEEN (
    SELECT ensembl_genomes_version FROM data_release WHERE is_current = 1
  ) - 1
  AND (
    SELECT ensembl_genomes_version FROM data_release WHERE is_current = 1
  )
  -- unique 'genebuild_method' values:
  -- ('full_genebuild', 'import', 'external_annotation_import', 'curated', 'anno', 'braker', 'projection_build')
  AND genome_annotation.value IN ('full_genebuild', 'anno', 'braker')
  -- result samples for development
  -- AND rand() <= 0.01
  -- AND rand() <= 0.001
;
