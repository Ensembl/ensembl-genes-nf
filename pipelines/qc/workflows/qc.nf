include { AGAT_METRICS } from '../subworkflows/agat/agat_stats.nf'
include { INTERPRO_SCAN } from '../subworkflows/interpro/interproscan.nf'
include { DIAMOND_PROTEIN_VALIDATION } from '../subworkflows/diamond/diamond_protein_validation.nf'
include { PREP_AA_FASTA } from '../modules/sequence/prep_aa_fasta.nf'

workflow QC {

    take:
        annotation_ch
        mode
        feature_levels_yaml
        database
        data_file_path
        diamond_reference_proteins
        diamond_reference_db

    main:
        // Each row supplies the shared annotation inputs. A precomputed protein
        // FASTA is used when present; otherwise it is derived with gffread.
        samples = annotation_ch.map { row ->
            def meta = [id: row[0], sample: row[0]]
            def protein = row.size() > 3 && row[3] ? file(row[3]) : null
            tuple(meta, file(row[1]), file(row[2]), protein)
        }

        gff_ch = samples.map { meta, _genome_fasta, gff3, _protein -> tuple(meta, gff3) }

        if (mode in ['agat', 'combined']) {
            AGAT_METRICS(
                gff_ch,
                feature_levels_yaml
            )
        }

        if (mode in ['interproscan', 'diamond', 'combined']) {
            supplied_proteins = samples
                .filter { _meta, _genome_fasta, _gff3, protein -> protein }
                .map { meta, _genome_fasta, _gff3, protein -> tuple(meta, protein) }

            samples_to_derive = samples
                .filter { _meta, _genome_fasta, _gff3, protein -> !protein }
                .map { meta, genome_fasta, gff3, _protein ->
                    tuple(meta, genome_fasta, gff3)
                }

            derived_proteins = PREP_AA_FASTA(samples_to_derive)
            protein_ch = supplied_proteins.mix(derived_proteins.aa_fasta)
        }

        if (mode in ['interproscan', 'combined']) {
            INTERPRO_SCAN(
                protein_ch,
                database,
                data_file_path
            )
        }

        if (mode in ['diamond', 'combined']) {
            DIAMOND_PROTEIN_VALIDATION(
                protein_ch,
                diamond_reference_proteins,
                diamond_reference_db
            )
        }
}
