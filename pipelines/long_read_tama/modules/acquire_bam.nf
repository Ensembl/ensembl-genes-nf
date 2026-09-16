process ACQUIRE_LONG_READ_BAM {
    tag "${meta.id}:${meta.classification}"
    label 'process_high_memory'

    input:
    tuple val(meta), val(expected_md5), val(expected_filename), val(source_uri), val(sidecar_uris), val(sidecar_md5s)
    val cache_dir

    output:
    tuple val(meta), path('downloaded/*'), emit: artifacts
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p downloaded "${cache_dir}/${meta.classification}/${expected_md5}"
    cache_file="${cache_dir}/${meta.classification}/${expected_md5}/${expected_filename}"
    lock="\${cache_file}.lock"
    while ! mkdir "\${lock}" 2>/dev/null; do sleep 10; done
    trap 'rmdir "\${lock}" 2>/dev/null || true' EXIT
    if [ ! -f "\${cache_file}" ]; then
      curl --fail --location --retry 3 --output "\${cache_file}.part" "${source_uri}"
      test "\$(md5sum "\${cache_file}.part" | awk '{print \$1}')" = "${expected_md5}"
      mv "\${cache_file}.part" "\${cache_file}"
    fi
    test "\$(md5sum "\${cache_file}" | awk '{print \$1}')" = "${expected_md5}" || { echo "BAM checksum mismatch for ${meta.id}" >&2; exit 1; }
    cp -p "\${cache_file}" "downloaded/${expected_filename}"
    IFS=';' read -ra URIS <<< "${sidecar_uris}"
    IFS=';' read -ra MD5S <<< "${sidecar_md5s}"
    if [ -n "${sidecar_uris}" ] && [ "\${#URIS[@]}" -ne "\${#MD5S[@]}" ]; then echo 'sidecar URI/checksum mismatch' >&2; exit 1; fi
    for i in "\${!URIS[@]}"; do
      side_name="\$(basename "\${URIS[i]}")"
      side_cache="${cache_dir}/${meta.classification}/\${MD5S[i]}/\${side_name}"
      mkdir -p "\$(dirname "\${side_cache}")"
      if [ ! -f "\${side_cache}" ]; then curl --fail --location --retry 3 --output "\${side_cache}.part" "\${URIS[i]}"; test "\$(md5sum "\${side_cache}.part" | awk '{print \$1}')" = "\${MD5S[i]}" || { echo "Sidecar checksum mismatch for ${meta.id}: \${side_name}" >&2; exit 1; }; mv "\${side_cache}.part" "\${side_cache}"; fi
      test "\$(md5sum "\${side_cache}" | awk '{print \$1}')" = "\${MD5S[i]}" || { echo "Cached sidecar checksum mismatch for ${meta.id}: \${side_name}" >&2; exit 1; }
      cp -p "\${side_cache}" downloaded/"\${side_name}"
    done
    printf 'run_accession\\texpected_md5\\tactual_md5\\tfilename\\n${meta.id}\\t${expected_md5}\\t${expected_md5}\\t${expected_filename}\\n' > downloaded/checksum.tsv
    printf '"%s":\\n    acquisition: curl-and-md5sum\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p downloaded
    touch downloaded/${expected_filename}
    touch downloaded/placeholder.pbi
    printf 'run_accession\\texpected_md5\\tactual_md5\\tfilename\\n${meta.id}\\tstub\\tstub\\t${expected_filename}\\n' > downloaded/checksum.tsv
    printf '"%s":\\n    acquisition: stub\\n' '${task.process}' > versions.yml
    """
}
