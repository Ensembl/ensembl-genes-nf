# MSA tools image

This image contains Python, curl, and the pinned CESAR2 `mafIndex` and
`mafExtract` binaries used by MSA reference setup and regional extraction.
Builds must be published by CI for linux/amd64, matching the Ensembl HPC
execution nodes, and production configuration must use the resulting manifest
digest. The upstream CESAR binaries are x86-64; a local tag is development-only.
