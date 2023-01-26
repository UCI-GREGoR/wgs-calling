rule create_exon_bedfile:
    """
    Extract exon annotations from gene annotation gtf
    and create annotated compressed bedfile for bcftools stats
    """
    input:
        "reference_data/references/{}/exons.gtf.gz".format(reference_build),
    output:
        bed_gz="results/bcftools_stats/exons.bed.gz",
        bed_gz_tbi="results/bcftools_stats/exons.bed.gz.tbi",
    benchmark:
        "results/performance_benchmarks/bcftools_stats/exons.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="1000",
        qname="small",
    shell:
        'gunzip -c {input} | awk -F"\\t" \'$3 == "exon" {{print $1"\\t"$4"\\t"$5}}\' | sort -k 1,1 -k 2,2g -k 3,3g | bgzip -c > {output.bed_gz} && '
        "tabix -s1 -b2 -e3 {output.bed_gz}"


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
        bed_gz="results/bcftools_stats/exons.bed.gz",
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
        "bcftools stats -E {input.bed_gz} -F {input.fasta} {input.vcf} > {output.stats}"
