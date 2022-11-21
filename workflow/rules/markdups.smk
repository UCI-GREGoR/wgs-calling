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
    conda:
        "../envs/bwamem2.yaml"
    threads: 4
    resources:
        h_vmem="16000",
        qname="small",
    shell:
        "samtools markdup --write-index --threads {threads} "
        "-f {output.score} -O BAM {output.bam}##idx##{output.bai}"
