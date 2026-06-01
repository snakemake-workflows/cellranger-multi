#!/usr/bin/env bash
set -euo pipefail

# set all the necessary conda paths and
# ensure they exist
CONDA_BIN="${CONDA_PREFIX}/bin"
mkdir -p ${CONDA_BIN}

CONDA_LIB="${CONDA_PREFIX}/lib"
mkdir -p ${CONDA_LIB}

MAIN_DIR=$( pwd )

# install cellranger
cd ${CONDA_LIB}
tar xzf ${CELLRANGER_TARBALL}
CELLRANGER_DIR=$( ls -d cellranger* )

ln -s ${CONDA_LIB}/${CELLRANGER_DIR}/cellranger ${CONDA_BIN}/cellranger

# disable telemetry
cellranger telemetry disable

# check that the cellranger executable is available and works
cellranger testrun --id=tiny --localmem=8

echo ${CONDA_BIN}/cellranger >${MAIN_DIR}/logs/cellranger/cellranger_conda_bin.txt