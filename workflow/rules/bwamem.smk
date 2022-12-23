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
        mem_mb="4000",
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
            suffix=bwa_ref_suffixes,
        ),
    params:
        exec_name=config["behaviors"]["aligner"],
    benchmark:
        "results/performance_benchmarks/bwa_index/{}/ref.fasta.tsv".format(reference_build)
    conda:
        "../envs/bwamem2.yaml" if (
        config["behaviors"]["aligner"] == "bwa-mem2"
        ) else "../envs/bwamem.yaml"
    threads: 1
    resources:
        mem_mb="64000",
        qname="small",
    shell:
        "{params.exec_name} index {input.fasta}"


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
            suffix=bwa_ref_suffixes + ["fai"],
        ),
    output:
        bam="results/aligned/{projectid}/{sampleid}_{lane}.bam",
        bai="results/aligned/{projectid}/{sampleid}_{lane}.bam.bai",
    benchmark:
        "results/performance_benchmarks/bwa_map_and_sort/{projectid}/{sampleid}_{lane}.tsv"
    params:
        exec_name=config["behaviors"]["aligner"],
        K="1000000",
        readgroup=lambda wildcards: "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:{}\\tPU:{}".format(
            "RG1", wildcards.sampleid, wildcards.sampleid, "Illumina", wildcards.sampleid
        ),
        tmpdir="temp",
    conda:
        "../envs/bwamem2.yaml" if (
        config["behaviors"]["aligner"] == "bwa-mem2"
        ) else "../envs/bwamem.yaml"
    threads: 12
    resources:
        mem_mb="500000",
        qname="large",
        tmpdir="temp",
    shell:
        "mkdir -p {params.tmpdir} && "
        '{params.exec_name} mem -t {threads} -Y -R "{params.readgroup}" -K {params.K} '
        "{input.bwa_fasta} {input.fastq1} {input.fastq2} | "
        "samtools sort -l 1 -m 4G -@ {threads} -T ${{TMPDIR}} -O BAM --write-index -o {output.bam}##idx##{output.bai}"
