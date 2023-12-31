rule bqsr_base_recalibrator:
    """
    Recalibrate base quality scores
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        index_files=expand(
            "reference_data/{aligner}/{genome}/ref.{suffix}",
            aligner=config["behaviors"]["aligner"],
            genome=reference_build,
            suffix=["dict", "fasta.fai"],
        ),
        known_indels="reference_data/bqsr/{}/ref.known.indels.vcf.gz".format(reference_build),
        known_indels_tbi="reference_data/bqsr/{}/ref.known.indels.vcf.gz.tbi".format(
            reference_build
        ),
        dbsnp138="reference_data/bqsr/{}/ref.dbsnp138.vcf".format(reference_build),
        dbsnp138_idx="reference_data/bqsr/{}/ref.dbsnp138.vcf.idx".format(reference_build),
    output:
        table="results/bqsr/{projectid}/{sampleid}.recal_table",
    params:
        tmpdir=tempDir,
        java_args=config_resources["gatk_bqsr_base_recalibrator"]["java_args"],
    benchmark:
        "results/performance_benchmarks/bqsr_base_recalibrator/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/gatk4.yaml" if not use_containers else None
    container:
        "{}/gatk4.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["gatk_bqsr_base_recalibrator"]["threads"]
    resources:
        mem_mb=config_resources["gatk_bqsr_base_recalibrator"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["gatk_bqsr_base_recalibrator"]["queue"], config_resources["queues"]
        ),
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" BaseRecalibrator '
        "--tmp-dir {params.tmpdir} "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-O {output.table} "
        "--known-sites {input.known_indels} "
        "--known-sites {input.dbsnp138}"


rule bqsr_apply_bqsr:
    """
    Apply BQSR to sample
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        index_files=expand(
            "reference_data/{aligner}/{genome}/ref.fasta.{suffix}",
            aligner=config["behaviors"]["aligner"],
            genome=reference_build,
            suffix=["dict", "fai"],
        ),
        table="results/bqsr/{projectid}/{sampleid}.recal_table",
    output:
        bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        bai="results/aligned_bams/{projectid}/{sampleid}.bai",
    params:
        tmpdir=tempDir,
        java_args=config_resources["gatk_bqsr_apply_bqsr"]["java_args"],
        use_bqsr=config["behaviors"]["bqsr"],
    benchmark:
        "results/performance_benchmarks/bqsr_apply_bqsr/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/gatk4.yaml" if not use_containers else None
    container:
        "{}/gatk4.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["gatk_bqsr_apply_bqsr"]["threads"]
    resources:
        mem_mb=config_resources["gatk_bqsr_apply_bqsr"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["gatk_bqsr_apply_bqsr"]["queue"], config_resources["queues"]
        ),
    shell:
        'if [[ "{params.use_bqsr}" == "False" ]] ; then cp {input.bam} {output.bam} && cp {input.bai} {output.bai} ; else '
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" ApplyBQSR '
        "--tmp-dir {params.tmpdir} "
        "-R {input.fasta} "
        "-I {input.bam} "
        "--bqsr-recal-file {input.table} "
        "-O {output.bam} ; "
        "fi"
