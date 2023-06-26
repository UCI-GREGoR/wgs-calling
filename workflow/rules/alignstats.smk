rule run_alignstats:
    """
    Run alignstats utility on post-markdups bams
    """
    input:
        bam="results/aligned_bams/{fileprefix}.bam",
        bai="results/aligned_bams/{fileprefix}.bai",
    output:
        json="results/alignstats/{fileprefix}.alignstats.json",
    benchmark:
        "results/performance_benchmarks/run_alignstats/{fileprefix}.tsv"
    params:
        min_qual=20,
    conda:
        "../envs/alignstats.yaml" if not use_containers else None
    container:
        "{}/alignstats.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["alignstats"]["threads"]
    resources:
        mem_mb=config_resources["alignstats"]["memory"],
        qname=rc.select_queue(config_resources["alignstats"]["queue"], config_resources["queues"]),
    shell:
        "alignstats -C -U "
        "-i {input.bam} "
        "-o {output.json} "
        "-j bam "
        "-v -P 10 -p 10 "
        "-q {params.min_qual}"


rule merge_alignstats:
    """
    Combine json-format alignstats output into a single big table
    """
    input:
        json=lambda wildcards: tc.construct_alignstats_targets(wildcards, manifest),
    output:
        yaml="results/alignstats/{projectid}/alignstats_summary_mqc.yaml",
    benchmark:
        "results/performance_benchmarks/merge_alignstats/{projectid}/alignstats_summary_mqc.tsv"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=rc.select_queue(config_resources["default"]["queue"], config_resources["queues"]),
    script:
        "../scripts/alignstats_json_to_yaml.py"
