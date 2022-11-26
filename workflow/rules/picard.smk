rule create_sequence_dictionary:
    """
    For a reference fasta, create a sequence dictionary (.dict extension)
    """
    input:
        "{prefix}fasta",
    output:
        "{prefix}fasta.dict",
    benchmark:
        "results/performance_benchmarks/create_sequence_dictionary/{prefix}fasta.tsv"
    params:
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        h_vmem="10000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
        'gatk --java-options "{params.java_args}" CreateSequenceDictionary '
        "-REFERENCE {input} "
        "-OUTPUT {output} "
        "--TMP_DIR {params.tmpdir}"


rule mark_duplicates:
    """
    Use samtools markdups to mark duplicates on aligned reads
    """
    input:
        bam="results/bwa-mem2/{fileprefix}.bwa2a.bam",
        bai="results/bwa-mem2/{fileprefix}.bwa2a.bam.bai",
    output:
        bam="results/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="results/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        score="results/markdups/{fileprefix}.mrkdup.score.txt",
    benchmark:
        "results/performance_benchmarks/mark_duplicates/{fileprefix}.tsv"
    params:
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        h_vmem="10000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
        'gatk --java-options "{params.java_args}" MarkDuplicates '
        "-INPUT {input.bam} "
        "-OUTPUT {output.bam} "
        "-METRICS_FILE {output.score} "
        "--CREATE_INDEX true "
        "--TMP_DIR {params.tmpdir} && "
        "mv results/markdups/{wildcards.fileprefix}.mrkdup.sort.bai {output.bai}"


rule picard_collectmultiplemetrics:
    """
    Run gatk version of picard CollectMultipleMetrics
    """
    input:
        bam="results/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="results/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
        dic="reference_data/references/{}/ref.fasta.dict".format(reference_build),
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
        stop_after="2000000",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        h_vmem="10000",
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
        "-STOP_AFTER {params.stop_after} "
        "-PROGRAM CollectAlignmentSummaryMetrics "
        "-PROGRAM CollectInsertSizeMetrics "
        "-PROGRAM QualityScoreDistribution "
        "-PROGRAM CollectSequencingArtifactMetrics "
        "-PROGRAM CollectQualityYieldMetrics "
        "-INCLUDE_UNPAIRED true "
        "-OUTPUT {params.outprefix}"


rule picard_collectgcbiasmetrics:
    """
    Run gatk version of picard CollectGcBiasMetrics
    """
    input:
        bam="results/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="results/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
        dic="reference_data/references/{}/ref.fasta.dict".format(reference_build),
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
        h_vmem="10000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
        'gatk --java-options "{params.java_args}" CollectGcBiasMetrics '
        "-INPUT {input.bam} "
        "-REFERENCE_SEQUENCE {input.fasta} "
        "-OUTPUT {output.metrics} "
        "-SUMMARY_OUTPUT {output.summary} "
        "-CHART_OUTPUT {output.pdf}"


rule picard_collectwgsmetrics:
    """
    Run gatk version of picard CollectWgsMetrics
    """
    input:
        bam="results/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="results/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
        dic="reference_data/references/{}/ref.fasta.dict".format(reference_build),
        intervals="reference_data/references/{}/ref.reportable-regions".format(reference_build),
    output:
        txt="results/collectwgsmetrics/{fileprefix}.picard.collect_wgs_metrics.txt",
    benchmark:
        "results/performance_benchmarks/picard_collectwgsmetrics/{fileprefix}.tsv"
    params:
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx6000m",
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        h_vmem="12000",
        qname="small",
        tmpdir="temp",
    shell:
        "mkdir -p temp/ && "
        'gatk --java-options "{params.java_args}" CollectWgsMetrics '
        "-INPUT {input.bam} "
        "-REFERENCE_SEQUENCE {input.fasta} "
        "-OUTPUT {output.txt} "
        "-INTERVALS {input.intervals}"
