#!/usr/bin/env bash
paste <(pigz -c -d $1 | paste - - - -) <(pigz -c -d ${1/_R1_/_R2_} | paste - - - -) | tr "\\t" "\\n"
