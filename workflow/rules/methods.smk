rule summarize_methods:
    """
    For end user consumption, generate a markdown description
    of the methods and user configuration settings in effect
    for an instance of the workflow.
    """
    input:
        manta_config=config["parameters"]["manta"]["config-ini"],
    output:
        markdown="results/reports/methods_summary.md",
    params:
        config=config,
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    script:
        "../scripts/summarize_workflow.py"
