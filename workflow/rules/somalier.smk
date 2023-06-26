rule somalier_extract:
    """
    Run somalier extract on a single bam.
    """
    input:
        bam="results/aligned_bams/{projectid}/{fileprefix}.bam",
        bai="results/aligned_bams/{projectid}/{fileprefix}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        sites_vcf="reference_data/somalier/{}/ref.sites.vcf.gz".format(reference_build),
    output:
        "results/somalier/{projectid}/extract/{fileprefix}.somalier",
    benchmark:
        "results/performance_benchmarks/somalier_extract/{projectid}/{fileprefix}.tsv"
    params:
        extract_dir="results/somalier/{projectid}/extract",
    conda:
        "../envs/somalier.yaml" if not use_containers else None
    container:
        "docker://brentp/somalier:v0.2.16" if use_containers else None
    threads: config_resources["somalier"]["threads"]
    resources:
        mem_mb=config_resources["somalier"]["memory"],
        qname=rc.select_queue(config_resources["somalier"]["queue"], config_resources["queues"]),
    shell:
        "somalier extract -d {params.extract_dir} "
        "--sites {input.sites_vcf} "
        "-f {input.fasta} "
        "{input.bam}"


rule somalier_relate:
    """
    Compute relatedness metrics on preprocessed alignment data with somalier.
    """
    input:
        somalier=lambda wildcards: tc.construct_somalier_extract_targets(wildcards, manifest),
        ped="results/somalier/{projectid}/relate/somalier.ped",
    output:
        html="results/somalier/{projectid}/relate/somalier.html",
        pairs="results/somalier/{projectid}/relate/somalier.pairs.tsv",
        samples="results/somalier/{projectid}/relate/somalier.samples.tsv",
    benchmark:
        "results/performance_benchmarks/somalier_relate/{projectid}.tsv"
    params:
        outprefix="results/somalier/{projectid}/relate/somalier",
    conda:
        "../envs/somalier.yaml" if not use_containers else None
    container:
        "docker://brentp/somalier:v0.2.16" if use_containers else None
    threads: config_resources["somalier"]["threads"]
    resources:
        mem_mb=config_resources["somalier"]["memory"],
        qname=rc.select_queue(config_resources["somalier"]["queue"], config_resources["queues"]),
    shell:
        "somalier relate --ped {input.ped} -o {params.outprefix} {input.somalier}"


rule somalier_build_pedfile:
    """
    Generate a pedfile for sex check for somalier.

    In the post 0.3.0 world, this is going to try to suck in self-reported sex
    information from the newly-generated sample ID linker. However, this information
    is only unreliably reported upstream, and as such this sexcheck will still be
    very incomplete. This is flagged to eventually be replaced with queries to retool,
    once that is actually implemented and operational.
    """
    input:
        linker="results/export/linker.tsv",
    output:
        ped="results/somalier/{projectid}/relate/somalier.ped",
        problems="results/somalier/{projectid}/discordant_annotations.tsv",
    benchmark:
        "results/performance_benchmarks/somalier_build_pedfile/{projectid}/somalier.tsv"
    params:
        subjectids=lambda wildcards: manifest.loc[
            manifest["projectid"] == wildcards.projectid, "sampleid"
        ].to_list(),
        ruid=lambda wildcards: wildcards.projectid,
        last_sample_sex=config["behaviors"]["assume-last-sample-sex"]
        if "assume-last-sample-sex" in config["behaviors"]
        else "unknown",
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=rc.select_queue(config_resources["default"]["queue"], config_resources["queues"]),
    script:
        "../scripts/construct_somalier_pedfile.py"
