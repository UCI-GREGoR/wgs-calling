rule samtools_index_fasta:
    """
    Given a fasta, generate its fai index file.
    """
    input:
        "{prefix}fasta",
    output:
        "{prefix}fasta.fai",
    benchmark:
        "results/performance_benchmarks/samtools_index_fasta/{prefix}fasta.fai.tsv"
    conda:
        "../envs/bwamem2.yaml"
    threads: 1
    resources:
        h_vmem="4000",
        qname="small",
    shell:
        "samtools faidx {input}"


rule bwa_index:
    """
    From a fasta file, run bwa index to generate annotation files
    """
    input:
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        index_files=expand(
            "reference_data/references/{genome}/ref.fasta.{suffix}",
            genome=reference_build,
            suffix=["amb", "ann", "bwt.2bit.64", "pac"],
        ),
    benchmark:
        "results/performance_benchmarks/bwa_index/{}/ref.fasta.tsv".format(reference_build)
    conda:
        "../envs/bwamem2.yaml"
    threads: 1
    resources:
        h_vmem="20000",
        qname="small",
    shell:
        "bwa index {input.fasta}"


rule bwa_map_and_sort:
    """
    Align fastqs to a reference genome
    """
    input:
        fastq1=lambda wildcards: tc.map_fastq_from_project_and_sample(
            wildcards, config, manifest, "R1"
        ),
        fastq2=lambda wildcards: tc.map_fastq_from_project_and_sample(
            wildcards, config, manifest, "R2"
        ),
        bwa_fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        bwa_other_files=expand(
            "reference_data/references/{genome}/ref.fasta.{suffix}",
            genome=reference_build,
            suffix=["ann", "amb", "bwt.2bit.64", "fai", "pac"],
        ),
    output:
        bam="results/bwa-mem2/{projectid}/{sampleid}.bwa2a.bam",
        bai="results/bwa-mem2/{projectid}/{sampleid}.bwa2a.bam.bai",
    benchmark:
        "results/performance_benchmarks/bwa_map_and_sort/{projectid}/{sampleid}.tsv"
    params:
        K="1000000",
        k="19",
        softclip_alts="",
        readgroup=lambda wildcards: "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:{}\\tPU:{}".format(
            "RG1", wildcards.sampleid, wildcards.sampleid, "Illumina", wildcards.sampleid
        ),
        tmpdir="temp",
    conda:
        "../envs/bwamem2.yaml"
    threads: 4
    resources:
        h_vmem="32000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p {params.tmpdir} && "
        'bwa-mem2 mem -t {threads} -K {params.K} -k {params.k} -Y -R "{params.readgroup}" {params.softclip_alts} '
        "{input.bwa_fasta} {input.fastq1} {input.fastq2} | "
        "samtools fixmate -@ {threads} -u -m - - | "
        "samtools sort -l 1 -m 2G -@ {threads} -T {params.tmpdir} -O BAM --write-index -o {output.bam}##idx##{output.bai}"
