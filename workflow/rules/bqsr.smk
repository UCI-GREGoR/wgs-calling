rule bqsr_base_recalibrator:
    """
    Recalibrate base quality scores
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        index_files=expand(
            "reference_data/references/{genome}/ref.fasta.{suffix}",
            genome=reference_build,
            suffix=["dict", "fai"],
        ),
        known_indels="reference_data/bqsr/{}/ref.known.indels.vcf.gz".format(reference_build),
        known_indels_tbi="reference_data/bqsr/{}/ref.known.indels.tbi".format(reference_build),
        dbsnp138="reference_data/bqsr/{}/ref.dbsnp138.vcf".format(reference_build),
        dbsnp138_idx="reference_data/bqsr/{}/ref.dbsnp138.vcf.idx".format(reference_build),
    output:
        table="results/bqsr/{projectid}/{sampleid}.recal_table",
    params:
        tmpdir="temp",
        java_args="-Xmx4000m -XX:+UseParallelGC -XX:ParallelGCThreads=16",
    benchmark:
        "results/performance_benchmarks/bqsr_base_recalibrator/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/gatk4.yaml"
    threads: 20
    resources:
        mem_mb="8000",
        qname="large",
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
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        index_files=expand(
            "reference_data/references/{genome}/ref.fasta.{suffix}",
            genome=reference_build,
            suffix=["dict", "fai"],
        ),
        table="results/bqsr/{projectid}/{sampleid}.recal_table",
    output:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
    params:
        tmpdir="temp",
        java_args="-Xmx10000m",
    benchmark:
        "results/performance_benchmarks/bqsr_base_recalibrator/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        mem_mb="20000",
        qname="large",
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" ApplyBQSR '
        "--tmp-dir {params.tmpdir} "
        "-R {input.fasta} "
        "-I {input.bam} "
        "--bqsr-recal-table {input.table} "
        "-O {output.bam}"
