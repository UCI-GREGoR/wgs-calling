rule merge_sv_vcfs:
    """
    Given SV callsets from various callers, combine into a single vcf
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sampleid}}.duphold-filtered.vcf.gz",
            toolname=config["behaviors"]["sv-callers"],
        ),
    output:
        "results/final/{projectid}/{sampleid}.sv.vcf.gz",
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
