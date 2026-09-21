#!/bin/bash

# ============================================================
# USER CONFIGURATION
# ============================================================

SAMPLE="p8_ee_Ztautau_ecm91"
FILENAME="events"
TAG="V3"

# Analysis script
ANALYSIS="analysisMakeTree.py"

# Fixed input base
INPUT_BASE="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/DelphesSampleWinter23/IDEA/"


# ============================================================
# Derived configuration
# ============================================================

INPUT_DIR="${INPUT_BASE}/${SAMPLE}"
OUTPUT_DIR="${PWD}/TREE_${SAMPLE}_${TAG}"
LOG_DIR="${OUTPUT_DIR}/condor_logs"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

echo "============================================================"
echo " Sample       : ${SAMPLE}"
echo " Filename     : ${FILENAME}"
echo " Input dir    : ${INPUT_DIR}"
echo " Output dir   : ${OUTPUT_DIR}"
echo " Analysis     : ${ANALYSIS}"
echo "============================================================"


# ============================================================
# Build input file list
# ============================================================

FILELIST="${OUTPUT_DIR}/input_files.txt"

rm -f "${FILELIST}"

for f in "${INPUT_DIR}"/"${FILENAME}"_*.root; do

    # Ignore the literal glob if there are no matching files
    [ -e "$f" ] || continue

    # Example:
    # events_123.root --> 123

    basename=$(basename "$f")

    number=${basename#${FILENAME}_}
    number=${number%.root}

    echo "${f} ${number}" >> "${FILELIST}"

done


NFILES=$(wc -l < "${FILELIST}")

echo "Found ${NFILES} files"

if [ "${NFILES}" -eq 0 ]; then
    echo "ERROR: no files found matching"
    echo "${INPUT_DIR}/${FILENAME}_*.root"
    exit 1
fi


# ============================================================
# Create Condor submit description
# ============================================================

SUBMIT_FILE="${OUTPUT_DIR}/condor.sub"

cat > "${SUBMIT_FILE}" <<EOF
universe = vanilla

executable = ${PWD}/run_condor_job_WIN23.sh

arguments = "\$(inputfile) \$(filenumber) ${OUTPUT_DIR} ${PWD}/${ANALYSIS}"

output = ${LOG_DIR}/job_\$(filenumber).out
error = ${LOG_DIR}/job_\$(filenumber).err
log = ${LOG_DIR}/job_\$(filenumber).log

+AccountingGroup = "group_u_CMST3.all"
+JobFlavour = "workday"

queue inputfile, filenumber from ${FILELIST}
EOF


# ============================================================
# Show submit configuration
# ============================================================

echo
echo "============================================================"
echo "Condor submit file:"
echo "${SUBMIT_FILE}"
echo "============================================================"
cat "${SUBMIT_FILE}"
echo "============================================================"
echo


# ============================================================
# Submit
# ============================================================

echo "Submitting ${NFILES} jobs to condorsc1.ciemat.es..."
echo

condor_submit -name condorsc1.ciemat.es "${SUBMIT_FILE}"
