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
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    script:
        "../scripts/summarize_workflow.py"
