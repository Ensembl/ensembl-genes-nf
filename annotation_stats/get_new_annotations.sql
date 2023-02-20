-- get new annotations


-- SELECT
--   *
-- FROM data_release
-- -- FROM data_release_database
-- -- WHERE rand() <= 0.2
-- -- WHERE rand() <= 0.1
-- -- WHERE rand() <= 0.01
-- ;


SELECT
  -- genome_annotation.genome_annotation_id AS 'genome_annotation_id',
  genome_annotation.type AS 'genome_annotation.type',
  genome_annotation.value AS 'genome_annotation.value',
  genome_database.genome_database_id AS 'genome_database_id',
  genome_database.dbname AS 'genome_database.dbname',
  genome.genome_id AS 'genome.genome_id',
  data_release.ensembl_version AS 'data_release.ensembl_version'
  -- COUNT(*)
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
  AND data_release.ensembl_version >= (
    SELECT ensembl_version FROM data_release WHERE is_current = 1
  ) - 1
  -- AND rand() <= 0.1
  -- AND rand() <= 0.01
  -- AND rand() <= 0.001
-- LIMIT 10
;
