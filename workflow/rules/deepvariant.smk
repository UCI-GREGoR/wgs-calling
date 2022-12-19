rule deepvariant_make_examples:
    """
    Run deepvariant make_examples in a hybrid
    embarrassingly parallel fashion.
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
        intervals="results/deepvariant/split_ranges/{splitnum}.ssv",
    output:
        temp(
            expand(
                "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config["parameters"]["deepvariant"]["number-shards"])
                ],
                shardmax=str(config["parameters"]["deepvariant"]["number-shards"]).rjust(5, "0"),
                suffix=["gz", "gz.example_info.json"],
            )
        ),
    benchmark:
        "results/performance_benchmarks/deepvariant_make_examples/{projectid}/{sampleid}.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config["parameters"]["deepvariant"]["number-shards"],
        ),
    container:
        "docker://google/deepvariant:{}".format(
            config["parameters"]["deepvariant"]["docker-version"]
        )
    threads: config["parameters"]["deepvariant"]["number-shards"]
    resources:
        mem_mb="{}".format(12000 * config["parameters"]["deepvariant"]["number-shards"]),
        qname="large",
        tmpdir="temp",
    shell:
        "mkdir -p temp && "
        "seq 0 $(({threads}-1)) | parallel -j{threads} --tmpdir ${{TMPDIR}} "
        "make_examples --mode calling "
        '--ref {input.fasta} --reads {input.bam} --regions "$(cat {input.intervals})" '
        "--examples {params.shard_string} --channels insert_size "
        "--task {{}}"


rule deepvariant_call_variants:
    """
    Run deepvariant call_variants in an
    embarrassingly parallel fashion.
    """
    input:
        expand(
            "results/deepvariant/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
            shardnum=[
                str(i).rjust(5, "0")
                for i in range(config["parameters"]["deepvariant"]["number-shards"])
            ],
            shardmax=str(config["parameters"]["deepvariant"]["number-shards"]).rjust(5, "0"),
            suffix=["gz", "gz.example_info.json"],
        ),
    output:
        gz=temp("results/deepvariant/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz"),
    benchmark:
        "results/performance_benchmarks/deepvariant_call_variants/{projectid}/{sampleid}.{splitnum}.tsv"
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
        mem_mb="{}".format(config["parameters"]["deepvariant"]["number-shards"] * 30000),
        qname="large",
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
        gz="results/deepvariant/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        vcf=temp(
            "results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.vcf.gz"
        ),
        tbi=temp(
            "results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.vcf.gz.tbi"
        ),
        html="results/deepvariant/{projectid}/postprocess_variants/{sampleid}.{splitnum}.visual_report.html",
    benchmark:
        "results/performance_benchmarks/deepvariant_postprocess_variants/{projectid}/{sampleid}.{splitnum}.tsv"
    container:
        "docker://google/deepvariant:{}".format(
            config["parameters"]["deepvariant"]["docker-version"]
        )
    threads: 1
    resources:
        mem_mb="32000",
        qname="large",
    shell:
        "postprocess_variants "
        "--ref {input.fasta} "
        "--infile {input.gz} "
        "--outfile {output.vcf}"


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
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "bcftools concat --threads {threads} -O z -o {output} {input}"
