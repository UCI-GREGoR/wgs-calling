## This file exists outside of the standard logic, at least for the time being.
## As such, let's try out some interesting things with the environment.

## For performance benchmarking, we'd like to know how many threads each rule
## has configured for its run. But that information can be hard to keep track of.
## We can *sort of* get it from the environment; but it introduces an order dependency
## on the inclusion of this file into the main workflow that is thoroughly antithetical
## to python. But I'm a C++ dev at heart, so eat it python.


def get_rule_threads(target_workflow):
    """
    Aggregate a dict of thread data for all loaded rules
    """
    rule_threads = {}
    for rule in target_workflow.rules:
        rule_threads[rule.name] = rule.resources["_cores"] if "_cores" in rule.resources else -1
    return rule_threads


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
            "deepvariant_make_examples",
            "deepvariant_call_variants",
            "deepvariant_postprocess_variants",
            "deepvariant_combine_regions",
            "duphold_apply",
            "duphold_run",
            "estimate_contamination",
            "index_vcf",
            "manta_configure",
            "manta_run",
            "manta_sort_output",
            "mark_duplicates",
            "merge_alignstats",
            "merge_sv_vcfs",
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
            "tiddit_sort_output",
        ],
        rule_threads=lambda wildcards: get_rule_threads(workflow),
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/aggregate_performance_metrics.Rmd"
