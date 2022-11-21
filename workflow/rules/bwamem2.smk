rule bwa_map_and_sort:
    """
    Align fastqs to a reference genome
    """
    input:
        fastq1=lambda wildcards: tc.map_fastq_from_project_and_sample(wildcards, manifest, "R1"),
        fastq2=lambda wildcards: tc.map_fastq_from_project_and_sample(wildcards, manifest, "R2"),
        bwa_fasta=config["bwa"]["grch37"]["bwa-fasta"],
    output:
        bam="results/bwa-mem2/{projectid}/{sampleid}.bwa2a.bam",
        bai="results/bwa-mem2/{projectid}/{sampleid}.bwa2a.bam.bai",
    params:
        K="60000000",
        k="19",
        softclip_alts="",
        readgroup=lambda wildcards: "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:\\tPU:{}".format(
            "RG1", wildcards.sampleid, wildcards.sampleid, "Illumina", wildcards.sampleid
        ),
        tmpdir="temp",
    conda:
        "../envs/bwamem2.yaml"
    threads: 12
    resources:
        h_vmem="24000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p {params.tmpdir} && "
        'bwa-mem2 mem -t {threads} -K {params.K} -k {params.k} -Y -R "{params.readgroup}" {params.softclip_alts} '
        "{input.bwa_fasta} {input.fastq1} {input.fastq2} | "
        "samtools fixmate -@ {threads} -u -m - - | "
        'samtools sort -l 1 -m 2G -@ {threads} -T {params.tmpdir} -O BAM --write-index -o "{output.bam}##{output.bai}"'
