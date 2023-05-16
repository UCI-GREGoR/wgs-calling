#!/usr/bin/env bash
if [[ "$#" -ne 2 ]] ; then
    echo "usage: ./export_data.bash {export-target-directory} {md5sum-results-file}"
    exit 1
fi
TARGET_DIR="$1"
RESULTS_FILE="$2"
RESULTS_FILE="$(readlink -m ${RESULTS_FILE})"

# determine jira ticket and flowcell from current path
# This is legacy support for the workflow when run on a particular
# infrastructure pre-April 2023. In the absence of the expected jira
# ticket or flowcell patterns, the workflow will still perform local
# data export, but the files will be emitted to the top level of
# the local export directory.
JIRA_TICKET=$(pwd | sed -r 's/.*(RT-[0-9]{4}).*/\1/')
FLOWCELL=$(pwd | sed -r 's/.*(RU[0-9]{5}).*/\1/')
TARGET_DIR="${TARGET_DIR}"
if [[ "${JIRA_TICKET}" == "RT-"* ]]  ; then
    TARGET_DIR="${TARGET_DIR}/${JIRA_TICKET}"
fi
if [[ "${FLOWCELL}" == "RU"* ]] ; then
    TARGET_DIR="${TARGET_DIR}/${FLOWCELL}"
fi
# create output target directory
mkdir -p "${TARGET_DIR}"
# copy targets
cp results/export/${FLOWCELL}/*vcf.gz* results/export/${FLOWCELL}/*cram* results/export/${FLOWCELL}/methods* "${TARGET_DIR}"
# return here later. this probably doesn't matter much
RETURN_LOCATION="$(pwd)"
cd "${TARGET_DIR}"
# run checksums
for file in $(ls *md5) ; do
    md5sum -c ${file} >> "${RESULTS_FILE}"
done
# go home
cd "${RETURN_LOCATION}"
