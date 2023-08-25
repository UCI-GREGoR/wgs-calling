#!/usr/bin/env bash
if [[ "$#" != "2" ]] ; then
    echo "usage: ./generate_manifest_upstream.bash s3://path/to/fastqs/ projectid"
    exit 1
fi
PROJECT_S3="$1"
PROJECT_ID="$2"
ALL_FASTQ=$(aws s3 ls --profile default "${PROJECT_S3}/" | awk -v s3path="${PROJECT_S3}" '{print s3path"/"$4}')
echo -e "projectid\tsampleid\tlane\tr1\tr2"
for filename in $(echo ${ALL_FASTQ} | sed 's/ /\n/g' | awk '/_R1.fastq.gz/') ; do
    SAMPLE_ID=$(echo ${filename} | awk -F"/" '{print $NF}' | cut -d"_" -f1)
    R1="${filename}"
    R2=$(echo ${filename} | sed 's/_R1.fastq.gz/_R2.fastq.gz/')
    aws s3 ls --profile default "${R2}" > /dev/null
    if [[ "$?" -eq "1" ]] ; then
	R2="NA"
    fi
    LANE="combined"
    if [[ "${SAMPLE_ID}" != "Undetermined" ]] ; then
	echo -e "${PROJECT_ID}\t${SAMPLE_ID}\t${LANE}\t${R1}\t${R2}" | sed 's|//|/|g' | sed 's|\ts3:/|\ts3://|g'
    fi
done
exit 0
