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
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        mem_mb="10000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
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
        bam=lambda wildcards: tc.get_bams_by_lane(wildcards, config, manifest, "bam"),
        bai=lambda wildcards: tc.get_bams_by_lane(wildcards, config, manifest, "bam.bai"),
    output:
        bam=temp("results/markdups/{projectid}/{sampleid}.mrkdup.bam"),
        score="results/markdups/{projectid}/{sampleid}.mrkdup.score.txt",
    benchmark:
        "results/performance_benchmarks/mark_duplicates/{projectid}/{sampleid}.tsv"
    params:
        tmpdir="temp",
        bamlist=lambda wildcards: " -INPUT ".join(
            tc.get_bams_by_lane(wildcards, config, manifest, "bam")
        ),
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx3000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 2
    resources:
        mem_mb="12000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
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
    benchmark:
        "results/performance_benchmarks/sort_bam/{prefix}.sort.bam.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    resources:
        mem_mb="8000",
        qname="small",
    shell:
        "samtools sort -@ {threads} -o {output.bam} -O bam {input.bam}"


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
        "../envs/samtools.yaml"
    threads: 4
    resources:
        mem_mb="8000",
        qname="small",
    shell:
        "samtools index -@ {threads} -b -o {output.bai} {input.bam}"


rule picard_collectmultiplemetrics:
    """
    Run gatk version of picard CollectMultipleMetrics
    """
    input:
        bam="results/bqsr/{fileprefix}.bam",
        bai="results/bqsr/{fileprefix}.bai",
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
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m",
        outprefix="results/collectmultiplemetrics/{fileprefix}.picard",
        extension=".txt",
        validation_stringency="LENIENT",
        metric_accumulation_level="SAMPLE",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        mem_mb="10000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
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
        bam="results/bqsr/{fileprefix}.bam",
        bai="results/bqsr/{fileprefix}.bai",
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
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        mem_mb="10000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
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
        bam="results/bqsr/{fileprefix}.bam",
        bai="results/bqsr/{fileprefix}.bai",
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
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx10000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        mem_mb="16000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
        'gatk --java-options "{params.java_args}" CollectWgsMetrics '
        "-INPUT {input.bam} "
        "-REFERENCE_SEQUENCE {input.fasta} "
        "-OUTPUT {output.txt} "
        "--TMP_DIR {params.tmpdir}"
