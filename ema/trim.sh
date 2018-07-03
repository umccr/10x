#!/bin/bash

FQ1=$1
FQ2=$2
FQ_OUT=$3
REPORT=$4
SAMPLE_ID=$5
THREADS=$6

# paired-> interleaved
atropos trim \
  --quality-base 33 --format fastq \
  -a 'A{200}' -a 'C{200}' -a 'G{200}' -a 'T{200}' \
  -A 'A{200}' -A 'C{200}' -A 'G{200}' -A 'T{200}' \
  --error-rate 0.3 \
  --overlap 8 \
  --no-default-adapters \
  --no-cache-adapters \
  -pe1 ${FQ1} \
  -pe2 ${FQ2} \
  --interleaved-output ${FQ_OUT} \
  --report-file ${REPORT} \
  --report-formats json \
  --sample-id ${SAMPLE_ID} \
  --quality-cutoff=5 \
  --minimum-length=25 \
  --threads ${THREADS}

# interleaved -> interleaved
#atropos trim \
#  --quality-base 33 --format fastq \
#  -a 'A{200}$' -a 'C{200}$' -a 'G{200}$' -a 'T{200}$' \
#  --overlap 8 --no-default-adapters --no-cache-adapters \
#  -A 'A{200}$' -A 'C{200}$' -A 'G{200}$' -A 'T{200}$' \
#  --interleaved-input ${FQ} \
#  --interleaved-output ${FQ_OUT} \
#  --report-file ${REPORT} \
#  --report-formats json \
#  --sample-id ${SAMPLE_ID} \
#  --quality-cutoff=5 \
#  --minimum-length=25
