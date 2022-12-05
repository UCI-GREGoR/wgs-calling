## This file exists outside of the standard logic, at least for the time being.
## As such, let's try out some interesting things with the environment.


def get_rule_threads(target_workflow):
    """
    Aggregate a dict of thread data for all loaded rules
    """
    rule_threads = {}
    for rule in target_workflow.rules:
        rule_threads[rule.name] = rule.resources["_cores"] if "_cores" in rule.resources else -1
    return rule_threads


def get_rule_names(target_workflow) -> list:
    """
    Aggregate a list of rule names for all loaded rules
    """
    return [rule.name for rule in target_workflow.rules]


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
        rules=lambda wildcards: get_rule_names(workflow),
        rule_threads=lambda wildcards: get_rule_threads(workflow),
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/aggregate_performance_metrics.Rmd"
