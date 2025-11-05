#!/usr/bin/env bash

( # keep everything logged
  any_missing=0
  echo "Checking these cellranger output files:"
  echo "${snakemake_output[*]}"
  for f in "${snakemake_output[@]}"
  do
    if [ ! -e "$f" ]
    then
      any_missing=1
      echo "Missing expected cellranger multi output: $f"
    fi
  done;
  if [ "$any_missing" -eq 1 ]
  then
    exit 1
  fi
  # if all exist, ensure timestamp of output is newer than the input
  touch "${snakemake_output[@]}"
) >${snakemake_log[0]} 2>&1