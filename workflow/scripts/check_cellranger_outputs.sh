#!/usr/bin/env bash

( # keep everything logged
  pwd
  any_missing=0
  echo "Will check files: ${snakemake_output[@]}"
  for f in "${snakemake_output[@]}"
  do
    echo "Now checking file: $f"
    ls -lh $f
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
  touch "${snakemake_output[@]}";
) >${snakemake_log[0]} 2>&1