rule create_sequence_dictionary:
    """
    For a reference fasta, create a sequence dictionary (.dict extension)

    Because GATK/bqsr is ridiculous, make two copies of the dict, so it's always available
    no matter how the downstream tool thinks it should be named.
    """
    input:
        "{prefix}fasta",
    output:
        standard="{prefix}fasta.dict",
        modified="{prefix}dict",
    benchmark:
        "results/performance_benchmarks/create_sequence_dictionary/{prefix}fasta.tsv"
    params:
        tmpdir=tempDir,
        java_args=config_resources["gatk_create_sequence_dictionary"]["java_args"],
    conda:
        "../envs/gatk4.yaml" if not use_containers else None
    container:
        "{}/gatk4.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["gatk_create_sequence_dictionary"]["threads"]
    resources:
        mem_mb=config_resources["gatk_create_sequence_dictionary"]["memory"],
        qname=rc.select_queue(
            config_resources["gatk_create_sequence_dictionary"]["queue"],
            config_resources["queues"],
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" CreateSequenceDictionary '
        "-REFERENCE {input} "
        "-OUTPUT {output.standard} "
        "--TMP_DIR {params.tmpdir} && "
        "cp {output.standard} {output.modified}"


rule mark_duplicates:
    """
    Use gatk/picard markdups to mark duplicates on aligned reads

    Note that the following tools are confirmed OK with marked but
    non-deduped reads:
    - deepvariant
    - manta
    - tiddit (? https://github.com/SciLifeLab/TIDDIT/blob/2bd94921ef2df4b08967e8f76fb10bc94730715d/tiddit/tiddit_signal.pyx#L160)
    - svaba

    Note that the following tools require deduping prior to use:
    - Sentieon DNAscope

    We are removing dups here, to save space and compute, and enable
    broader tool selection.  Keep an eye on this moving forward
    in case we ever want to change this around.
    """
    input:
        bam=lambda wildcards: tc.get_bams_by_lane(wildcards, checkpoints, config, manifest, "bam"),
        bai=lambda wildcards: tc.get_bams_by_lane(
            wildcards, checkpoints, config, manifest, "bam.bai"
        ),
    output:
        bam=temp("results/markdups/{projectid}/{sampleid}.mrkdup.bam"),
        score="results/markdups/{projectid}/{sampleid}.mrkdup.score.txt",
    benchmark:
        "results/performance_benchmarks/mark_duplicates/{projectid}/{sampleid}.tsv"
    params:
        tmpdir=tempDir,
        bamlist=lambda wildcards: " -INPUT ".join(
            tc.get_bams_by_lane(wildcards, checkpoints, config, manifest, "bam")
        ),
        java_args=config_resources["gatk_mark_duplicates"]["java_args"],
    conda:
        "../envs/gatk4.yaml" if not use_containers else None
    container:
        "{}/gatk4.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["gatk_mark_duplicates"]["threads"]
    resources:
        mem_mb=config_resources["gatk_mark_duplicates"]["memory"],
        qname=rc.select_queue(
            config_resources["gatk_mark_duplicates"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" MarkDuplicates '
        "-INPUT {params.bamlist} "
        "-OUTPUT {output.bam} "
        "-REMOVE_DUPLICATES true "
        "-METRICS_FILE {output.score} "
        "--CREATE_INDEX false "
        "--TMP_DIR {params.tmpdir}"


rule sort_bam:
    input:
        bam="{prefix}.bam",
    output:
        bam="{prefix}.sort.bam",
    params:
        tmpdir=tempDir,
    benchmark:
        "results/performance_benchmarks/sort_bam/{prefix}.sort.bam.tsv"
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/bwa.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools_sort"]["threads"]
    resources:
        mem_mb=config_resources["samtools_sort"]["memory"],
        qname=rc.select_queue(
            config_resources["samtools_sort"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        "samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} -O bam {input.bam}"


rule samtools_create_bai:
    """
    From a sorted bam file, create a bai-format index
    """
    input:
        bam="results/{prefix}.sort.bam",
    output:
        bai="results/{prefix}.sort.bam.bai",
    benchmark:
        "results/performance_benchmarks/samtools_create_bai/{prefix}.sort.tsv"
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/bwa.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools"]["threads"]
    resources:
        mem_mb=config_resources["samtools"]["memory"],
        qname=rc.select_queue(config_resources["samtools"]["queue"], config_resources["queues"]),
    shell:
        "samtools index -@ {threads} -b -o {output.bai} {input.bam}"


rule picard_collectmultiplemetrics:
    """
    Run gatk version of picard CollectMultipleMetrics
    """
    input:
        bam="results/aligned_bams/{fileprefix}.bam",
        bai="results/aligned_bams/{fileprefix}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        dic="reference_data/{}/{}/ref.fasta.dict".format(
            config["behaviors"]["aligner"], reference_build
        ),
    output:
        expand(
            "results/collectmultiplemetrics/{{fileprefix}}.picard.{suffix}",
            suffix=[
                "alignment_summary_metrics.txt",
                "base_distribution_by_cycle_metrics.txt",
                "bait_bias_summary_metrics.txt",
                "error_summary_metrics.txt",
                "insert_size_metrics.txt",
                "pre_adapter_summary_metrics.txt",
                "quality_by_cycle_metrics.txt",
                "quality_distribution_metrics.txt",
                "quality_yield_metrics.txt",
            ],
        ),
    benchmark:
        "results/performance_benchmarks/picard_collectmultiplemetrics/{fileprefix}.tsv"
    params:
        tmpdir=tempDir,
        java_args=config_resources["gatk_collectmultiplemetrics"]["java_args"],
        outprefix="results/collectmultiplemetrics/{fileprefix}.picard",
        extension=".txt",
        validation_stringency="LENIENT",
        metric_accumulation_level="SAMPLE",
    conda:
        "../envs/gatk4.yaml" if not use_containers else None
    container:
        "{}/gatk4.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["gatk_collectmultiplemetrics"]["threads"]
    resources:
        mem_mb=config_resources["gatk_collectmultiplemetrics"]["memory"],
        qname=rc.select_queue(
            config_resources["gatk_collectmultiplemetrics"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" CollectMultipleMetrics '
        "-INPUT {input.bam} "
        "-REFERENCE_SEQUENCE {input.fasta} "
        "-FILE_EXTENSION {params.extension} "
        "-VALIDATION_STRINGENCY {params.validation_stringency} "
        "-METRIC_ACCUMULATION_LEVEL {params.metric_accumulation_level} "
        "-LEVEL {params.metric_accumulation_level} "
        "-PROGRAM CollectAlignmentSummaryMetrics "
        "-PROGRAM CollectInsertSizeMetrics "
        "-PROGRAM QualityScoreDistribution "
        "-PROGRAM CollectSequencingArtifactMetrics "
        "-PROGRAM CollectQualityYieldMetrics "
        "-INCLUDE_UNPAIRED true "
        "-OUTPUT {params.outprefix} "
        "--TMP_DIR {params.tmpdir}"


rule picard_collectgcbiasmetrics:
    """
    Run gatk version of picard CollectGcBiasMetrics
    """
    input:
        bam="results/aligned_bams/{fileprefix}.bam",
        bai="results/aligned_bams/{fileprefix}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        dic="reference_data/{}/{}/ref.fasta.dict".format(
            config["behaviors"]["aligner"], reference_build
        ),
    output:
        metrics="results/collectgcbiasmetrics/{fileprefix}.picard.gc_bias_metrics.txt",
        summary="results/collectgcbiasmetrics/{fileprefix}.picard.gc_bias_metrics_summary.txt",
        pdf="results/collectgcbiasmetrics/{fileprefix}.picard.gc_bias_metrics_chart.pdf",
    benchmark:
        "results/performance_benchmarks/picard_collectgcbiasmetrics/{fileprefix}.tsv"
    params:
        tmpdir=tempDir,
        java_args=config_resources["gatk_collectgcbiasmetrics"]["java_args"],
    conda:
        "../envs/gatk4.yaml" if not use_containers else None
    container:
        "{}/gatk4.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["gatk_collectgcbiasmetrics"]["threads"]
    resources:
        mem_mb=config_resources["gatk_collectgcbiasmetrics"]["memory"],
        qname=rc.select_queue(
            config_resources["gatk_collectgcbiasmetrics"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" CollectGcBiasMetrics '
        "-INPUT {input.bam} "
        "-REFERENCE_SEQUENCE {input.fasta} "
        "-OUTPUT {output.metrics} "
        "-SUMMARY_OUTPUT {output.summary} "
        "-CHART_OUTPUT {output.pdf} "
        "--TMP_DIR {params.tmpdir}"


rule picard_collectwgsmetrics:
    """
    Run gatk version of picard CollectWgsMetrics
    """
    input:
        bam="results/aligned_bams/{fileprefix}.bam",
        bai="results/aligned_bams/{fileprefix}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        dic="reference_data/{}/{}/ref.fasta.dict".format(
            config["behaviors"]["aligner"], reference_build
        ),
    output:
        txt="results/collectwgsmetrics/{fileprefix}.picard.collect_wgs_metrics.txt",
    benchmark:
        "results/performance_benchmarks/picard_collectwgsmetrics/{fileprefix}.tsv"
    params:
        tmpdir=tempDir,
        java_args=config_resources["gatk_collectwgsmetrics"]["java_args"],
    conda:
        "../envs/gatk4.yaml" if not use_containers else None
    container:
        "{}/gatk4.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["gatk_collectwgsmetrics"]["threads"]
    resources:
        mem_mb=config_resources["gatk_collectwgsmetrics"]["memory"],
        qname=rc.select_queue(
            config_resources["gatk_collectwgsmetrics"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" CollectWgsMetrics '
        "-INPUT {input.bam} "
        "-REFERENCE_SEQUENCE {input.fasta} "
        "-OUTPUT {output.txt} "
        "--TMP_DIR {params.tmpdir}"
