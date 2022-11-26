#!/usr/bin/env bash
if [[ "$#" != "2" ]] ; then
    echo "usage: ./generate_manifest.bash /path/to/project/directory projectid"
    exit 1
fi
PROJECT_DIR="$1"
PROJECT_ID="$2"
ALL_FASTQ=$(ls "${PROJECT_DIR}"/output/*fastq.gz)
echo -e "projectid\tsampleid\tlane\tr1\tr2"
for filename in $(echo ${ALL_FASTQ} | sed 's/ /\n/g' | awk '/_R1_/') ; do
    SAMPLE_ID=$(echo ${filename} | sed 's|/|\t|g' | awk -F'\t' '{print $NF}' | sed 's/_/\n/' | awk 'NR == 1')
    R1="${filename}"
    R2=$(echo ${filename} | sed 's/_R1_/_R2_/')
    if [[ ! -f ${R2} ]] ; then
	R2="NA"
    fi
    LANE=$(echo ${filename} | sed 's|/|\t|g' | awk -F'\t' '{print $NF}' | sed 's/_/\t/g' | awk -F'\t' '{print $2}')
    if [[ "${SAMPLE_ID}" != "Undetermined" ]] ; then
	echo -e "${PROJECT_ID}\t${SAMPLE_ID}\t${LANE}\t${R1}\t${R2}" | sed 's|//|/|g'
    fi
done
exit 0
