#!/bin/bash

# ============================================================
# CONFIGURATION
# ============================================================

SAMPLE="p8_ee_Ztautau_ecm91"
FILENAME="events"
#"out_sim_edm4hep"

INPUT_BASE="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC//DelphesSampleWinter23/IDEA/"

INPUT_DIR="${INPUT_BASE}/${SAMPLE}"

echo $INPUT_DIR

# Files you want to run interactively
#FILES="1"

#FILES="115197259 000143148 000149335 000302188 000436491 000705878 000818249 000858512 001048812 001422495 001513079 001977576 002087779 002145500 002310115 002474453 002600645 002912243 003195252 003208683 003296487 003604878 003625215 003673558 004161082 004167492 004278484 004393106 004419858 004495983 004697112 004839082 004987034 005314140 005438434 005713622 006072859 006364260 006561104 006667445"

# Maximum number of files to process

FILES=""

for f in "$INPUT_DIR"/"${FILENAME}"_*.root; do
    [ -e "$f" ] || continue

    n="${f##*/}"           # events_120124599.root
    n="${n#${FILENAME}_}"  # 120124599.root
    n="${n%.root}"         # 120124599

    FILES="$FILES $n"
done

FILES=$(echo "$FILES" | xargs)

echo "FILES=\"$FILES\""

# Output tag
TAG="V2"

# ============================================================

OUTPUT_DIR="TREE_${SAMPLE}_${TAG}"
mkdir -p "${OUTPUT_DIR}"

# Setup Key4HEP
#source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08

echo "Sample:     ${SAMPLE}"
echo "Input:      ${INPUT_DIR}"
echo "Output:     ${OUTPUT_DIR}"
echo "Files:      ${FILES}"
echo

for N in ${FILES}; do

    INPUT="${INPUT_DIR}/${FILENAME}_${N}.root"
    OUTPUT="${OUTPUT_DIR}/tree_${N}.root"

    echo "Checking existing output: $OUTPUT"
    ls -l "$OUTPUT" 2>/dev/null

    if [[ -e "$OUTPUT" ]]; then
        echo "Output already exists, skipping: $OUTPUT"
        continue
    fi

    echo "============================================================"
    echo "Processing file ${N}"
    echo "Input : ${INPUT}"
    echo "Output: ${OUTPUT}"
    echo "============================================================"

    if [ ! -f "${INPUT}" ]; then
        echo "WARNING: ${INPUT} does not exist. Skipping."
        continue
    fi

    python analysisMakeTree.py \
        -f "${INPUT}" \
        -o "${OUTPUT}" \
        -hl "False" -w "True"


    if [ $? -ne 0 ]; then
        echo "ERROR processing file ${N}"
    else
        echo "Finished file ${N}"
    fi

done

echo
echo "Done."
