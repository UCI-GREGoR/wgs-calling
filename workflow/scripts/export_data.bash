#!/usr/bin/env bash
if [[ "$#" -ne 2 ]] ; then
    echo "usage: ./export_data.bash {export-target-directory} {md5sum-results-file}"
    exit 1
fi
TARGET_DIR="$1"
RESULTS_FILE="$2"
RESULTS_FILE="$(readlink -m ${RESULTS_FILE})"

# determine jira ticket and flowcell from current path
JIRA_TICKET=$(pwd | sed -r 's/.*(RT-[0-9]{4}).*/\1/')
FLOWCELL=$(pwd | sed -r 's/.*(RU[0-9]{5}).*/\1/')
if [[ "${JIRA_TICKET}" != "RT-"* ]] || [[ "${FLOWCELL}" != "RU"* ]] ; then
    echo "unable to extract jira ticket and flowcell identifier from pwd $(pwd)"
    exit 2
fi
# create output target directory
TARGET_DIR="${TARGET_DIR}/${JIRA_TICKET}/${FLOWCELL}"
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
