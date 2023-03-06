rule bqsr_base_recalibrator:
    """
    Recalibrate base quality scores, but with the spark implementation and run locally.

    Not sure exactly how this is gonna go but I think it's worth a shot.
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
        tmpdir="temp",
        java_opts="-XX:+UseParallelGC -XX:ParallelGCThreads=2",
        extra="--known-sites reference_data/bqsr/{}/ref.known.indels.vcf.gz".format(reference_build),
    benchmark:
        "results/performance_benchmarks/bqsr_base_recalibrator/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/gatk4.yaml"
    threads: 8
    resources:
        mem_mb="16000",
        qname="large",
    wrapper:
        "v1.22.0/bio/gatk/baserecalibratorspark"


rule bqsr_apply_bqsr:
    """
    Apply BQSR to sample, but with the spark implementation and run locally.

    Not sure exactly how this is gonna go but I think it's worth a shot.
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
        recal_table="results/bqsr/{projectid}/{sampleid}.recal_table",
    output:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
    params:
        tmpdir="temp",
    benchmark:
        "results/performance_benchmarks/bqsr_apply_bqsr_spark/{projectid}/{sampleid}.tsv"
    threads: 8
    resources:
        mem_mb="20000",
        qname="large",
    wrapper:
        "v1.22.0/bio/gatk/applybqsrspark"
