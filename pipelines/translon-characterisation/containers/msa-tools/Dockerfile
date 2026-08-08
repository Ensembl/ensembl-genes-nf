ARG CESAR2_COMMIT=ef1828c0cb783a1a3c26d5750cdd5682428d2ab3
FROM alpine/git:2.49.1 AS source
ARG CESAR2_COMMIT
RUN git clone https://github.com/hillerlab/CESAR2.0.git /src \
    && git -C /src checkout "${CESAR2_COMMIT}" \
    && rm -rf /src/.git

FROM quay.io/biocontainers/ucsc-bedtobigbed@sha256:910c3c3ad7fa0e8ac8930fda4b7b0fb4a47f58339752478d4443c64e2440ad40 AS ucsc

FROM python:3.12.8-slim-bookworm@sha256:8859bd6ca943079262c27e38b7119cdacede77c463139a15651dd340087a6cc9
COPY --from=ucsc /usr/local/bin/bedToBigBed /usr/local/bin/bedToBigBed
COPY --from=ucsc /usr/local/lib/ /usr/local/lib/
COPY --from=source /src/tools/mafIndex /usr/local/bin/mafIndex
COPY --from=source /src/tools/mafExtract /usr/local/bin/mafExtract
COPY bedSort /usr/local/bin/bedSort
COPY --from=source /src/LICENSE /usr/local/share/licenses/CESAR2.0.LICENSE
RUN chmod 0755 /usr/local/bin/mafIndex /usr/local/bin/mafExtract /usr/local/bin/bedSort

LABEL org.opencontainers.image.source="https://github.com/hillerlab/CESAR2.0"
LABEL org.opencontainers.image.revision="ef1828c0cb783a1a3c26d5750cdd5682428d2ab3"
