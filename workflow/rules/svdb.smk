rule merge_sv_vcfs:
    """
    Given SV callsets from various callers, combine into a single vcf
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sampleid}}.{toolname}.duphold-filtered.vcf.gz",
            toolname=config["behaviors"]["sv-callers"],
        ),
    output:
        temp("results/final/{projectid}/{sampleid}.sv.svdb-raw.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/merge_sv_vcfs/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/svdb.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
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
        "results/final/{projectid}/{sampleid}.sv.svdb-raw.vcf.gz",
    output:
        "results/final/{projectid}/{sampleid}.sv.vcf.gz",
    params:
        bcftools_filter="|".join(
            [".*" for i in range(config["behaviors"]["sv-ensemble-min-count"])]
        ),
    benchmark:
        "results/performance_benchmarks/ensemble_sv_vcf/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "bcftools filter -i \"INFO/svdb_origin ~ '{params.bcftools_filter}' \" -O z -o {output} {input}"
