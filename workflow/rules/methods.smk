rule summarize_methods:
    """
    For end user consumption, generate an html description
    of the methods and user configuration settings in effect
    for an instance of the workflow.
    """
    input:
        jinja_template="workflow/scripts/summarize_workflow.jinja",
        manta_config=config["parameters"]["manta"]["config-ini"],
    output:
        html="{prefix}/methods_summary.html",
    params:
        config=config,
    conda:
        "../envs/python_jinja.yaml"
    container:
        "{}/python_jinja.sif".format(apptainer_images)
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=rc.select_queue(config_resources["default"]["queue"], config_resources["queues"]),
    script:
        "../scripts/summarize_workflow.py"
