rule merge_sv_vcfs:
    """
    Given SV callsets from various callers, combine into a single vcf
    """
    input:
        lambda wildcards: expand(
            "results/{toolname}/{{projectid}}/{{sampleid}}.{toolname}.duphold-filtered.vcf.gz",
            toolname=config["behaviors"]["sv-endpoints"][wildcards.endpoint]["sv-callers"],
        ),
    output:
        temp("results/final/{projectid}/{sampleid}.sv.{endpoint}.svdb-raw.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/merge_sv_vcfs/{projectid}/{sampleid}.{endpoint}.tsv"
    conda:
        "../envs/svdb.yaml" if not use_containers else None
    container:
        "{}/svdb.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["svdb"]["threads"]
    resources:
        mem_mb=config_resources["svdb"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["svdb"]["queue"], config_resources["queues"]
        ),
    shell:
        "svdb --merge --vcf {input} | "
        "sed 's/ID=PL,Number=G,Type=Integer/ID=PL,Number=G,Type=Float/' | "
        "sed 's/ID=GQ,Number=1,Type=Integer/ID=GQ,Number=1,Type=String/' | "
        "bgzip -c > {output}"


rule ensemble_sv_vcf:
    """
    Given an svdb-merged version of SV calls, select variants
    that are flagged as present in at least some user-configurable
    number of SV caller inputs
    """
    input:
        "results/final/{projectid}/{sampleid}.sv.{endpoint}.svdb-raw.vcf.gz",
    output:
        "results/final/{projectid}/{sampleid}.sv.{endpoint}.vcf.gz",
    params:
        bcftools_filter_count=lambda wildcards: "INFO/FOUNDBY >= {}".format(
            config["behaviors"]["sv-endpoints"][wildcards.endpoint]["sv-ensemble"]["min-count"]
        ),
        bcftools_filter_sources=" & INFO/svdb_origin ~ '"
        + "' & INFO/svdb_origin ~ '".join(
            config["behaviors"]["sv-endpoints"][wildcards.endpoint]["sv-ensemble"][
                "required-callers"
            ]
        )
        + "'"
        if "required-callers"
        in config["behaviors"]["sv-endpoints"][wildcards.endpoint]["sv-ensemble"]
        else "",
    benchmark:
        "results/performance_benchmarks/ensemble_sv_vcf/{projectid}/{sampleid}.{endpoint}.tsv"
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
        'bcftools filter -i "{params.bcftools_filter_count} {params.bcftools_filter_sources}" -O z -o {output} {input}'


rule summarize_sv_variant_sources:
    """
    For a raw merge output from svdb, use bcftools to extract reasonable summary information about
    the variant from its info content, with the goal of getting this information into an Rmd report
    """
    input:
        "results/final/{projectid}/{sampleid}.sv.{endpoint}.svdb-raw.vcf.gz",
    output:
        "results/reports/sv_data/{projectid}/{sampleid}.sv.{endpoint}.svdb-raw.tsv",
    benchmark:
        "results/performance_benchmarks/summarize_sv_variant_sources/{projectid}/{sampleid}.{endpoint}.tsv"
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
        "bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/svdb_origin\\n' {input} > {output}"
