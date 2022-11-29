rule deepvariant_make_examples:
    """
    Run deepvariant make_examples in an
    embarrassingly parallel fashion.
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
        intervals="results/deepvariant/split_ranges/{splitnum}.ssv",
    output:
        expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
            shardnum=[i + 1 for i in range(config["parameters"]["deepvariant"]["number-shards"])],
            shardmax=config["parameters"]["deepvariant"]["number-shards"],
            suffix=["gz", "gz.run_info.pbtxt"],
        ),
    params:
        shard_string=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config["parameters"]["deepvariant"]["number-shards"],
        ),
    container:
        "docker://google/deepvariant:{}".format(
            config["parameters"]["deepvariant"]["docker-version"]
        )
    threads: 1
    resources:
        h_vmem="2000",
        qname="small",
    shell:
        "make_examples --mode calling "
        '--ref {input.fasta} --reads {input.bam} --regions "$(cat {input.intervals})" '
        "--examples {params.shard_string} --channels insert_size "
        "--task {wildcards.shardnum}"


rule deepvariant_call_variants:
    """
    Run deepvariant call_variants in an
    embarrassingly parallel fashion.
    """
    input:
        expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
            shardnum=[i for i in range(config["parameters"]["deepvariant"]["number-shards"])],
            shardmax=config["parameters"]["deepvariant"]["number-shards"],
            suffix=["gz", "gz.run_info.pbtxt"],
        ),
    output:
        gz="results/deepvariant/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz",
    params:
        shard_string=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config["parameters"]["deepvariant"]["number-shards"],
        ),
        docker_model="/opt/models/wgs/model.ckpt",
    container:
        "docker://google/deepvariant:{}".format(
            config["parameters"]["deepvariant"]["docker-version"]
        )
    threads: config["parameters"]["deepvariant"]["number-shards"]
    resources:
        h_vmem="{}".format(config["parameters"]["deepvariant"]["number-shards"] * 4000),
        qname="small",
    shell:
        "call_variants "
        "--outfile {output.gz} "
        "--examples {params.shard_string} "
        '--checkpoint "{params.docker_model}"'


rule deepvariant_postprocess_variants:
    """
    Run deepvariant postprocess_variants in an
    embarrassingly parallel fashion.
    """
    input:
        "results/deepvariant/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz",
    output:
        "results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.vcf.gz",
    container:
        "docker://google/deepvariant:{}".format(
            config["parameters"]["deepvariant"]["docker-version"]
        )
    threads: 1
    resources:
        h_vmem="32000",
        qname="small",
    shell:
        "postprocess_variants "
        "--ref {input.fasta} "
        "--infile {input} "
        "--outfile {output}"


rule deepvariant_combine_regions:
    """
    Combine per-region deepvariant vcfs
    into a single mega vcf.
    """
    input:
        expand(
            "results/deepvariant/{{projectid}}/postprocess_variants/{{sampleid}}.{splitnum}.vcf.gz",
            splitnum=[i + 1 for i in range(octopus_num_intervals)],
        ),
    output:
        "results/deepvariant/{projectid}/{sampleid}.sorted.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        h_vmem="4000",
        qname="small",
    shell:
        "bcftools concat --threads {threads} -O z -o {output} {input}"
