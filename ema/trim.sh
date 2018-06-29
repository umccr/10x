#!/bin/bash

FQ=$1
FQ_OUT=$2
REPORT=$3
SAMPLE_ID=$4

atropos trim \
  --quality-base 33 --format fastq \
  -a 'A{200}$' -a 'C{200}$' -a 'G{200}$' -a 'T{200}$' \
  --overlap 8 --no-default-adapters --no-cache-adapters \
  -A 'A{200}$' -A 'C{200}$' -A 'G{200}$' -A 'T{200}$' \
  -l ${FQ} \
  -o >(bgzip --threads 1 -c > ${FQ_OUT}) \
  --report-file ${REPORT} \
  --report-formats json \
  --sample-id ${SAMPLE_ID} \
  --quality-cutoff=5 \
  --minimum-length=25 \
  --nextseq-trim=25
