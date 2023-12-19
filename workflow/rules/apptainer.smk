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
            config["parameters"]["deepvariant"]["docker-version"]
        ),
    params:
        image_version=config["parameters"]["deepvariant"]["docker-version"],
        tmpdir=tempDir,
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["apptainer"]["threads"]
    resources:
        mem_mb=config_resources["apptainer"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["apptainer"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        "APPTAINER_TMPDIR=$(realpath {params.tmpdir}) apptainer build --disable-cache {output} docker://google/deepvariant:{params.image_version}"
