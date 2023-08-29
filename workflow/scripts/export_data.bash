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
TARGET_DIR="${TARGET_DIR}"
if [[ "${JIRA_TICKET}" == "RT-"* ]]  ; then
    TARGET_DIR="${TARGET_DIR}/${JIRA_TICKET}"
fi
# create output target directory and copy targets
for projectid in $(find results/export -mindepth 1 -maxdepth 1 -type d -print | sed 's|results/export/||') ; do
    mkdir -p ${TARGET_DIR}/${projectid}
    find results/export/${projectid} \( -name "*vcf.gz*" -o -name "*cram*" -o -name "*crai*" -o -name "methods*" -o -name "*zip" \) -exec cp {} "${TARGET_DIR}/${projectid}" \;
    RETURN_LOCATION="$(pwd)"
    cd "${TARGET_DIR}/${projectid}"
    # run checksums
    for file in $(find . -name "*md5" -print) ; do
        md5sum -c ${file} >> "${RESULTS_FILE}"
    done
    # go home
    cd "${RETURN_LOCATION}"
done
# handle the situation where no files are exported
touch "${RESULTS_FILE}"
