rule mark_duplicates:
    """
    Use samtools markdups to mark duplicates on aligned reads
    """
    input:
        bam="{pathprefix}/bwa-mem2/{fileprefix}.bwa2a.bam",
        bai="{pathprefix}/bwa-mem2/{fileprefix}.bwa2a.bam.bai",
    output:
        bam="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        score="{pathprefix}/markdups/{fileprefix}.mrkdup.score.txt",
    params:
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        h_vmem="10000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
        'gatk --java-options "{params.java_args}" MarkDuplicates '
        "-INPUT {input.bam} "
        "-OUTPUT {output.bam} "
        "-METRICS_FILE {output.score} "
        "--CREATE_INDEX true "
        "--TMP_DIR {params.tmpdir} && "
        "mv {wildcards.pathprefix}/markdups/{wildcards.fileprefix}.mrkdup.sort.bai {output.bai}"
