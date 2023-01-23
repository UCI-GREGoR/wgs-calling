rule bcftools_stats:
    """
    Run bcftools stats on a tabix-indexed vcf file

    This rule is patterned in such a way that it is expected
    to be run against pre- and post-filter vcfs
    """
    input:
        vcf="results/{prefix}.vcf.gz",
        tbi="results/{prefix}.vcf.gz.tbi",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        stats="results/{prefix}.vcf.stats",
    benchmark:
        "results/performance_benchmarks/bcftools_stats/{prefix}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="8000",
        qname="small",
    shell:
        "bcftools stats -F {input.fasta} {input.vcf} > {output.stats}"
