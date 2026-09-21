#!/bin/bash

set -e

INPUT="$1"
JOBID="$2"
OUTPUT_DIR="$3"
CHANNEL="$4"
PMAX="$5"
MMAX="$6"

FINAL_OUTPUT="${OUTPUT_DIR}/histos_${JOBID}.root"

echo "============================================================"
echo "Host         : $(hostname)"
echo "Input        : ${INPUT}"
echo "Job ID       : ${JOBID}"
echo "Channel      : ${CHANNEL}"
echo "Final output : ${FINAL_OUTPUT}"
echo "Condor scratch: ${_CONDOR_SCRATCH_DIR}"
echo "============================================================"

# ------------------------------------------------------------------
# If this input has already been processed, don't redo it.
# Useful for resubmissions.
# ------------------------------------------------------------------

if [ -f "${FINAL_OUTPUT}" ]; then
    echo "Output already exists, skipping:"
    echo "  ${FINAL_OUTPUT}"
    exit 0
fi

# ------------------------------------------------------------------
# Work ONLY in the private Condor sandbox.
# makeLFVHistos_selection.py and reconstructZMass.py have been
# transferred here by Condor.
# ------------------------------------------------------------------

cd "${_CONDOR_SCRATCH_DIR}"

echo "Working directory:"
pwd

echo "Transferred files:"
ls -lh

source /cvmfs/sw.hsf.org/key4hep/setup.sh

echo "Python:"
which python
python --version

# Local output exists only inside this job's private scratch area.
LOCAL_OUTPUT="histos_${JOBID}.root"

python analysisGENFromTREE.py \
    -f "${INPUT}" \
    -o "${LOCAL_OUTPUT}" 

# ------------------------------------------------------------------
# Safety checks before copying
# ------------------------------------------------------------------

if [ ! -s "${LOCAL_OUTPUT}" ]; then
    echo "ERROR: output file was not produced or is empty"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Copy to a temporary name in the destination filesystem first.
# mv is then atomic, so another process never sees a half-written ROOT file.
TMP_OUTPUT="${OUTPUT_DIR}/.histos_${JOBID}.root.$$.tmp"

cp "${LOCAL_OUTPUT}" "${TMP_OUTPUT}"
mv "${TMP_OUTPUT}" "${FINAL_OUTPUT}"

echo "============================================================"
echo "Finished successfully:"
ls -lh "${FINAL_OUTPUT}"
echo "============================================================"
