#!/bin/bash

set -e

# ============================================================
# USER CONFIGURATION
# ============================================================

BASE="/nfs/cms/cepeda/FCC/AFBPlots/makePlots"

SAMPLE="p8_ee_Ztautau_ecm91_V3"
TAG="Win23_1498_V2"

INPUT_DIR="${BASE}/../makeTrees/TREE_${SAMPLE}"

OUTPUT_DIR="${BASE}/HIST_${SAMPLE}_${TAG}"

SCHEDD="condorsc1.ciemat.es"

# ============================================================

mkdir -p "${OUTPUT_DIR}"
mkdir -p condor_logs

rm -f joblist.txt

N=0

for INPUT in "${INPUT_DIR}"/*.root; do

    [ -e "${INPUT}" ] || continue

    # e.g. tree_157.root -> 157
    BASENAME=$(basename "${INPUT}" .root)
    JOBID="${BASENAME#tree_}"

    echo "${INPUT} ${JOBID} ${OUTPUT_DIR} "  >> joblist.txt

    N=$((N+1))
done

echo "============================================================"
echo "Sample      : ${SAMPLE}"
echo "Input dir   : ${INPUT_DIR}"
echo "Output dir  : ${OUTPUT_DIR}"
echo "Jobs        : ${N}"
echo "============================================================"

if [ "${N}" -eq 0 ]; then
    echo "ERROR: no ROOT files found"
    exit 1
fi

condor_submit -name "${SCHEDD}" submit_makeAFBHistos.sub
