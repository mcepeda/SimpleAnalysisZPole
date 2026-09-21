#!/bin/bash

INPUTFILE=$1
FILENUMBER=$2
OUTPUTDIR=$3
ANALYSIS=$4

echo "============================================================"
echo "Starting job"
echo "Host       : $(hostname)"
echo "Input      : ${INPUTFILE}"
echo "File number: ${FILENUMBER}"
echo "Output dir : ${OUTPUTDIR}"
echo "Analysis   : ${ANALYSIS}"
echo "============================================================"

# ------------------------------------------------------------
# Key4HEP environment
# ------------------------------------------------------------

source /cvmfs/sw.hsf.org/spackages6/key4hep-stack/2022-12-23/x86_64-centos7-gcc11.2.0-opt/ll3gi/setup.sh

echo
echo "Python:"
which python
python --version
echo

# Condor scratch directory
cd "${_CONDOR_SCRATCH_DIR}"

OUTPUTFILE="tree_${FILENUMBER}.root"

echo "Running:"
echo "python ${ANALYSIS} -f ${INPUTFILE} -o ${OUTPUTFILE} --hl F"
echo

python "${ANALYSIS}" \
    -f "${INPUTFILE}" \
    -o "${OUTPUTFILE}" \
    -hl "False" -w "True"

STATUS=$?

if [ ${STATUS} -ne 0 ]; then
    echo "ERROR: analysis failed with status ${STATUS}"
    exit ${STATUS}
fi

if [ ! -f "${OUTPUTFILE}" ]; then
    echo "ERROR: expected output ${OUTPUTFILE} was not created"
    exit 1
fi

# ------------------------------------------------------------
# Copy output
# ------------------------------------------------------------

echo
echo "Copying ${OUTPUTFILE} -> ${OUTPUTDIR}/"

cp "${OUTPUTFILE}" "${OUTPUTDIR}/"

STATUS=$?

if [ ${STATUS} -ne 0 ]; then
    echo "ERROR: output copy failed"
    exit ${STATUS}
fi

echo
echo "Successfully finished:"
echo "${OUTPUTDIR}/${OUTPUTFILE}"
