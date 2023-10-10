rule apptainer_deepvariant:
    """
    There have been sporadic instances of snakemake not being
    able to run deepvariant's official container, due to
    snakemake's reliance on an explicit `bash` command within
    the container. The workaround for this is an obnoxious
    direct call to apptainer, but which requires a sif file directly.
    """
    output:
        "results/apptainer_images/deepvariant_{}.sif".format(
            config["deepvariant"]["docker-version"]
        ),
    params:
        outdir="results/apptainer_images",
        image_version=config["deepvariant"]["docker-version"],
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "apptainer pull --dir {params.outdir} docker://google/deepvariant:{params.image_version}"
