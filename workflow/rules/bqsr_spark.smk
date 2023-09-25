rule bqsr_base_recalibrator_spark:
    """
    Recalibrate base quality scores, but with the spark implementation and run locally.

    This is *apparently* functional, but the spark implementation itself is flagged as
    in permanent beta and not ready for production.
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        ref="reference_data/{}/{}/ref.fasta".format(config["behaviors"]["aligner"], reference_build),
        dict=expand(
            "reference_data/{aligner}/{genome}/ref.dict",
            aligner=config["behaviors"]["aligner"],
            genome=reference_build,
        ),
        fai=expand(
            "reference_data/{aligner}/{genome}/ref.fasta.fai",
            aligner=config["behaviors"]["aligner"],
            genome=reference_build,
        ),
        extra="reference_data/bqsr/{}/ref.known.indels.vcf.gz".format(reference_build),
        extra_tbi="reference_data/bqsr/{}/ref.known.indels.vcf.gz.tbi".format(reference_build),
        known="reference_data/bqsr/{}/ref.dbsnp138.vcf".format(reference_build),
        known_idx="reference_data/bqsr/{}/ref.dbsnp138.vcf.idx".format(reference_build),
    output:
        recal_table="results/bqsr/{projectid}/{sampleid}.recal_table",
    params:
        tmpdir=tempDir,
        java_opts=config_resources["gatk_bqsr_base_recalibrator"]["java_args"],
        extra="--known-sites reference_data/bqsr/{}/ref.known.indels.vcf.gz".format(reference_build),
    benchmark:
        "results/performance_benchmarks/bqsr_base_recalibrator_spark/{projectid}/{sampleid}.tsv"
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
    wrapper:
        "v2.6.0/bio/gatk/baserecalibratorspark"


if config["behaviors"]["bqsr"]:

    rule bqsr_apply_bqsr_spark:
        """
        Apply BQSR to sample, but with the spark implementation and run locally.

        This is not functional, and the spark implementation itself is flagged as
        in permanent beta and not ready for production.
        """
        input:
            bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
            bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
            ref="reference_data/{}/{}/ref.fasta".format(
                config["behaviors"]["aligner"], reference_build
            ),
            dict=expand(
                "reference_data/{aligner}/{genome}/ref.dict",
                aligner=config["behaviors"]["aligner"],
                genome=reference_build,
            ),
            fai=expand(
                "reference_data/{aligner}/{genome}/ref.fasta.fai",
                aligner=config["behaviors"]["aligner"],
                genome=reference_build,
            ),
            recal_table="results/bqsr/{projectid}/{sampleid}.recal_table",
        output:
            bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        params:
            tmpdir=tempDir,
            java_args=config_resources["gatk_bqsr_apply_bqsr"]["java_args"],
            use_bqsr=config["behaviors"]["bqsr"],
        benchmark:
            "results/performance_benchmarks/bqsr_apply_bqsr_spark/{projectid}/{sampleid}.tsv"
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
        wrapper:
            "v2.6.0/bio/gatk/applybqsrspark"

else:

    rule copy_markdups_output:
        """
        Respond to user request to disable bqsr by circumventing bqsr entirely
        """
        input:
            bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
            bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        output:
            bam="results/aligned_bams/{projectid}/{sampleid}.bam",
            bai="results/aligned_bams/{projectid}/{sampleid}.bai",
        shell:
            "cp {input.bam} {output.bam} && "
            "cp {input.bai} {output.bai}"
