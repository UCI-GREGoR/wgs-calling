rule performance_benchmarks:
    """
    Combine performance benchmark data for various rules
    into a single report.
    """
    input:
        r_resources="workflow/scripts/aggregate_performance_metrics.R",
    output:
        "results/performance_benchmarks/performance_benchmarks.html",
    params:
        parent_dir="results/performance_benchmarks",
        rules=[
            "bwa_map_and_sort",
            "copy_fastqs",
            "estimate_contamination",
            "manta_configure",
            "manta_run",
            "mark_duplicates",
            "merge_alignstats",
            "picard_collectgcbiasmetrics",
            "picard_collectmultiplemetrics",
            "picard_collectwgsmetrics",
            "run_alignstats",
            "run_fastp",
            "run_mosdepth",
            "run_multiqc_alignment",
            "run_multiqc_fastq",
            "somalier_extract",
            "somalier_relate",
            "tiddit_run",
        ],
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/aggregate_performance_metrics.Rmd"
