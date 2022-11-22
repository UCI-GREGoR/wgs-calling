rule mark_duplicates:
    """
    Use samtools markdups to mark duplicates on aligned reads
    """
    input:
        bam="{pathprefix}/bwa-mem2/{fileprefix}.bwa2a.bam",
        bai="{pathprefix}/bwa-mem2/{fileprefix}.bwa2a.bam.bai",
    output:
        bam="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        score="{pathprefix}/markdups/{fileprefix}.mrkdup.score.txt",
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
        "mv {wildcards.pathprefix}/markdups/{wildcards.fileprefix}.mrkdup.sort.bai {output.bai}"


rule picard_collectmultiplemetrics:
    """
    Run gatk version of picard CollectMultipleMetrics
    """
    input:
        bam="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        fasta=config["references"]["grch37"]["fasta"],
        fai=config["references"]["grch37"]["fasta-index"],
        dic=config["references"]["grch37"]["fasta-dict"],
    output:
        expand(
            "{{pathprefix}}/collectmultiplemetrics/{{fileprefix}}.picard.{suffix}",
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
    params:
        tmpdir="temp",
        java_args="-Djava.io.tmpdir=temp/ -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m",
        outprefix="{pathprefix}/collectmultiplemetrics/{fileprefix}.picard",
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
        bam="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        fasta=config["references"]["grch37"]["fasta"],
        fai=config["references"]["grch37"]["fasta-index"],
        dic=config["references"]["grch37"]["fasta-dict"],
    output:
        metrics="{pathprefix}/collectgcbiasmetrics/{fileprefix}.picard.gc_bias_metrics.txt",
        summary="{pathprefix}/collectgcbiasmetrics/{fileprefix}.picard.gc_bias_metrics_summary.txt",
        pdf="{pathprefix}/collectgcbiasmetrics/{fileprefix}.picard.gc_bias_metrics_chart.pdf",
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
        "-CHART_OUTPUT {output.chart}"
