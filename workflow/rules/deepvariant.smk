localrules:
    caller_scatter_tasks,


rule caller_scatter_tasks:
    """
    For scatter-gather: take input set of range annotations
    and prepare for individual caller runs
    """
    input:
        lambda wildcards: config[wildcards.caller][reference_build]["calling-ranges"],
    output:
        "results/{caller}/split_ranges/{splitnum}.ssv",
    benchmark:
        "results/performance_benchmarks/{caller}/split_ranges/{splitnum}.tsv"
    shell:
        "cat $(awk 'NR == {wildcards.splitnum}' {input}) | tr '\\n' ' ' > {output}"


rule deepvariant_make_examples:
    """
    Run deepvariant make_examples in a hybrid
    embarrassingly parallel fashion.
    """
    input:
        bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        bai="results/aligned_bams/{projectid}/{sampleid}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        intervals="results/deepvariant/split_ranges/{splitnum}.ssv",
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["parameters"]["deepvariant"]["docker-version"]
        ),
    output:
        temp(
            expand(
                "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
                shardnum=[
                    str(i).rjust(5, "0") for i in range(config_resources["deepvariant"]["threads"])
                ],
                shardmax=str(config_resources["deepvariant"]["threads"]).rjust(5, "0"),
                suffix=["gz", "gz.example_info.json"],
            )
        ),
        temp(
            expand(
                "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
                shardnum=[
                    str(i).rjust(5, "0") for i in range(config_resources["deepvariant"]["threads"])
                ],
                shardmax=str(config_resources["deepvariant"]["threads"]).rjust(5, "0"),
            )
        ),
    benchmark:
        "results/performance_benchmarks/deepvariant_make_examples/{projectid}/{sampleid}.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deepvariant"]["threads"],
        ),
        gvcf_string=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deepvariant"]["threads"],
        ),
        tmpdir="/tmp",
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["deepvariant"]["threads"]
    resources:
        mem_mb=config_resources["deepvariant"]["make_examples_memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["deepvariant"]["queue"], config_resources["queues"]
        ),
        tmpdir="/tmp",
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "mkdir -p {params.tmpdir} && '
        "seq 0 $(({threads}-1)) | parallel -j{threads} --tmpdir {params.tmpdir} "
        "make_examples --mode calling "
        "--ref {input.fasta} "
        "--reads {input.bam} "
        '--regions "$(cat {input.intervals})" '
        "--examples {params.shard_string} --channels insert_size "
        "--gvcf {params.gvcf_string} "
        '--task {{}}"'


rule deepvariant_call_variants:
    """
    Run deepvariant call_variants in an
    embarrassingly parallel fashion.
    """
    input:
        infiles=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
            shardnum=[
                str(i).rjust(5, "0") for i in range(config_resources["deepvariant"]["threads"])
            ],
            shardmax=str(config_resources["deepvariant"]["threads"]).rjust(5, "0"),
            suffix=["gz", "gz.example_info.json"],
        ),
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["parameters"]["deepvariant"]["docker-version"]
        ),
    output:
        gz=temp("results/deepvariant/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz"),
    benchmark:
        "results/performance_benchmarks/deepvariant_call_variants/{projectid}/{sampleid}.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deepvariant"]["threads"],
        ),
        docker_model="/opt/models/wgs/model.ckpt",
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["deepvariant"]["threads"]
    resources:
        mem_mb=config_resources["deepvariant"]["call_variants_memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["deepvariant"]["queue"], config_resources["queues"]
        ),
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "'
        "call_variants "
        "--outfile {output.gz} "
        "--examples {params.shard_string} "
        '--checkpoint \\"{params.docker_model}\\""'


rule deepvariant_postprocess_variants:
    """
    Run deepvariant postprocess_variants in an
    embarrassingly parallel fashion.
    """
    input:
        gz="results/deepvariant/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        gvcf=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
            shardnum=[
                str(i).rjust(5, "0") for i in range(config_resources["deepvariant"]["threads"])
            ],
            shardmax=str(config_resources["deepvariant"]["threads"]).rjust(5, "0"),
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["parameters"]["deepvariant"]["docker-version"]
        ),
    output:
        vcf=temp(
            "results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.vcf.gz"
        ),
        gvcf=temp(
            "results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.g.vcf.gz"
        ),
        tbi=temp(
            "results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.vcf.gz.tbi"
        ),
        html="results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.visual_report.html",
    params:
        gvcf_string=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deepvariant"]["threads"],
        ),
    benchmark:
        "results/performance_benchmarks/deepvariant_postprocess_variants/{projectid}/{sampleid}.{splitnum}.tsv"
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: 1
    resources:
        mem_mb=config_resources["deepvariant"]["postprocess_variants_memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["deepvariant"]["queue"], config_resources["queues"]
        ),
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "'
        "postprocess_variants "
        "--ref {input.fasta} "
        "--infile {input.gz} "
        "--nonvariant_site_tfrecord_path {params.gvcf_string} "
        "--gvcf_outfile {output.gvcf} "
        '--outfile {output.vcf}"'


rule deepvariant_combine_regions:
    """
    Combine per-region deepvariant vcfs
    into a single mega vcf.
    """
    input:
        expand(
            "results/deepvariant/{{projectid}}/postprocess_variants/{{sampleid}}.{splitnum}.vcf.gz",
            splitnum=[i + 1 for i in range(caller_num_intervals)],
        ),
    output:
        "results/deepvariant/{projectid}/{sampleid}.sorted.vcf.gz",
    benchmark:
        "results/performance_benchmarks/deepvariant_combine_regions/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bcftools concat --threads {threads} -O z -o {output} {input}"


use rule deepvariant_combine_regions as deepvariant_combine_gvcfs with:
    input:
        expand(
            "results/deepvariant/{{projectid}}/postprocess_variants/{{sampleid}}.{splitnum}.g.vcf.gz",
            splitnum=[i + 1 for i in range(caller_num_intervals)],
        ),
    output:
        "results/deepvariant/{projectid}/{sampleid}.sorted.g.vcf.gz",
    benchmark:
        "results/performance_benchmarks/deepvariant_combine_gvcfs/{projectid}/{sampleid}.tsv"
