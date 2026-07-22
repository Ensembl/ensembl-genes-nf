#!/usr/bin/env nextflow
/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

nextflow.enable.dsl=2

include { GENMAP_INDEX } from 'nf-core/genmap/index'
include { GENMAP_MAP } from 'nf-core/genmap/map'
include { UCSC_BEDGRAPHTOBIGWIG } from 'nf-core/ucsc/bedgraphtobigwig'

def ensembl_species_name(String species) {
    def tokens = species.split('_')
    ([tokens[0].capitalize()] + tokens.drop(1).collect { token -> token.toLowerCase() }).join('_')
}

def ensembl_reference_url(String species, String release) {
    def release_path = release?.toLowerCase() == 'current' ? 'current' : "release-${release}"
    "https://ftp.ensembl.org/pub/${release_path}/fasta/${species}/dna/"
}

process FETCH_REFERENCE {
    tag "${meta.id}:download"
    label 'fetch_file'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    tuple val(meta), val(reference_url)

    output:
    tuple val(meta), path('reference.fa'), emit: reference
    path 'reference_source.txt', emit: source

    script:
    """
    reference_url='${reference_url}'
    if [[ "\$reference_url" == */ ]]; then
        reference_name=\$(curl --fail --location "\$reference_url" \\
            | sed -n 's/.*href="\\([^"]*\\.dna\\.toplevel\\.fa\\.gz\\)".*/\\1/p' \\
            | head -n 1)
        test -n "\$reference_name"
        reference_url="\${reference_url}\${reference_name}"
    fi
    curl --fail --location --retry 3 --retry-delay 5 \\
        --output reference.download "\$reference_url"
    if gzip -t reference.download 2>/dev/null; then
        gzip --decompress --stdout reference.download > reference.fa
    else
        mv reference.download reference.fa
    fi
    test -s reference.fa
    printf '%s\\n' "\$reference_url" > reference_source.txt
    """

    stub:
    """
    touch reference.fa reference_source.txt
    """
}

process USE_LOCAL_REFERENCE {
    tag "${meta.id}:local-reference"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path('reference.fa'), emit: reference
    path 'reference_source.txt', emit: source

    script:
    """
    if gzip -t '${reference}' 2>/dev/null; then
        gzip --decompress --stdout '${reference}' > reference.fa
    else
        cp '${reference}' reference.fa
    fi
    test -s reference.fa
    printf '%s\\n' 'local:${reference}' > reference_source.txt
    """

    stub:
    """
    touch reference.fa reference_source.txt
    """
}

process MAKE_REGIONS {
    tag "${meta.id}:regions"
    label 'default'
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: 'reference.sizes'

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path('regions.bed'), emit: regions
    path 'reference.sizes', emit: sizes

    script:
    """
    : > regions.bed
    : > reference.sizes
    awk '
        /^>/ {
            if (name) {
                print name "\\t0\\t" seq_length >> "regions.bed"
                print name "\\t" seq_length >> "reference.sizes"
            }
            name = substr(\$1, 2)
            seq_length = 0
            next
        }
        { seq_length += length(\$0) }
        END {
            if (name) {
                print name "\\t0\\t" seq_length >> "regions.bed"
                print name "\\t" seq_length >> "reference.sizes"
            }
        }
    ' '${reference}'
    test -s regions.bed
    test -s reference.sizes
    """

    stub:
    """
    touch regions.bed
    touch reference.sizes
    """
}

workflow {
    def species = params.species ?: 'reference'

    if (!params.outdir) error 'Missing required parameter: --outdir'
    if (!params.kmer.toString().isInteger() || params.kmer.toInteger() <= 0) {
        error '--kmer must be a positive integer'
    }
    if (!params.mismatches.toString().isInteger() || params.mismatches.toInteger() < 0) {
        error '--mismatches must be a non-negative integer'
    }
    if (params.reference_url && params.reference) {
        error 'Use only one of --reference_url and --reference'
    }
    if (!params.reference_url && !params.reference && !params.species) {
        error 'Provide --reference_url, --reference, or --species'
    }

    def source_id
    if (params.species) {
        source_id = "${params.species.toString().toLowerCase()}.release${params.release}"
    } else if (params.reference_url) {
        source_id = params.reference_url.toString().tokenize('/').last()
    } else {
        source_id = file(params.reference, checkIfExists: true).name
    }
    source_id = source_id.replaceAll(/(?i)\.(fa|fasta|fna)(\.gz)?$/, '')
    source_id = source_id.replaceAll('[^A-Za-z0-9_.-]', '_')
    def output_id = "${source_id}.k${params.kmer}.e${params.mismatches}"
    def meta = [id: output_id, species: species, source: source_id]

    def reference
    if (params.reference_url) {
        def reference_input = channel.value(tuple(meta, params.reference_url.toString()))
        FETCH_REFERENCE(reference_input)
        reference = FETCH_REFERENCE.out.reference
    } else if (params.reference) {
        def reference_input = channel.value(tuple(meta, file(params.reference, checkIfExists: true)))
        USE_LOCAL_REFERENCE(reference_input)
        reference = USE_LOCAL_REFERENCE.out.reference
    } else {
        def url = ensembl_reference_url(params.species.toString().toLowerCase(), params.release.toString())
        def reference_input = channel.value(tuple(meta, url))
        FETCH_REFERENCE(reference_input)
        reference = FETCH_REFERENCE.out.reference
    }

    def indexed = GENMAP_INDEX(reference)
    def regions = MAKE_REGIONS(reference)
    def mappability = GENMAP_MAP(
        indexed.index,
        regions.regions
    )
    UCSC_BEDGRAPHTOBIGWIG(
        mappability.bedgraph,
        regions.sizes
    )
}
