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
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=rc.select_queue(config_resources["default"]["queue"], config_resources["queues"]),
    shell:
        "samtools faidx {input}"


rule bwa_index:
    """
    From a fasta file, run bwa index to generate annotation files
    """
    input:
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
    output:
        index_files=expand(
            "reference_data/{aligner}/{genome}/ref.fasta.{suffix}",
            aligner=config["behaviors"]["aligner"],
            genome=reference_build,
            suffix=aligner_index_suffixes[config["behaviors"]["aligner"]],
        ),
    params:
        exec_name=config["behaviors"]["aligner"],
    benchmark:
        "results/performance_benchmarks/{}_index/{}/ref.fasta.tsv".format(
            config["behaviors"]["aligner"], reference_build
        )
    conda:
        "../envs/{}.yaml".format(config["behaviors"]["aligner"])
    threads: config_resources["bwa_index"]["threads"]
    resources:
        mem_mb=config_resources["bwa_index"]["memory"],
        qname=rc.select_queue(config_resources["bwa_index"]["queue"], config_resources["queues"]),
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
        bwa_fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        bwa_other_files=expand(
            "reference_data/{aligner}/{genome}/ref.fasta.{suffix}",
            aligner=config["behaviors"]["aligner"],
            genome=reference_build,
            suffix=aligner_index_suffixes[config["behaviors"]["aligner"]] + ["fai"],
        ),
    output:
        bam=temp("results/aligned/{projectid}/{sampleid}_{lane}.bam"),
        bai=temp("results/aligned/{projectid}/{sampleid}_{lane}.bam.bai"),
    benchmark:
        "results/performance_benchmarks/bwa_map_and_sort/{projectid}/{sampleid}_{lane}.tsv"
    params:
        exec_name=config["behaviors"]["aligner"],
        K=config["parameters"][config["behaviors"]["aligner"]]["K"],
        readgroup=lambda wildcards: "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:{}\\tPU:{}.{}.{}".format(
            "RG" + wildcards.lane,
            wildcards.sampleid,
            wildcards.sampleid,
            "ILLUMINA",
            wildcards.projectid,
            wildcards.lane,
            wildcards.sampleid,
        ),
        tmpdir=tempDir,
    conda:
        lambda wildcards: "../envs/{}.yaml".format(config["behaviors"]["aligner"])
    threads: config_resources["bwa_map_and_sort"]["threads"]
    resources:
        mem_mb=config_resources["bwa_map_and_sort"]["memory"],
        qname=rc.select_queue(
            config_resources["bwa_map_and_sort"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        '{params.exec_name} mem -t {threads} -Y -R "{params.readgroup}" -K {params.K} '
        "{input.bwa_fasta} {input.fastq1} {input.fastq2} | "
        "samtools sort -l 1 -m 4G -@ {threads} -T ${{TMPDIR}} -O BAM --write-index -o {output.bam}##idx##{output.bai}"
