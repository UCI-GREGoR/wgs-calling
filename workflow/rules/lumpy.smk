rule lumpy_run:
    """
    Run lumpy, using the smoove interface as recommended by the lumpy maintainers
    """
    input:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        bed="reference_data/lumpy/{}/ref.exclude.bed".format(reference_build),
    output:
        vcf="results/lumpy/{projectid}/{sampleid}.lumpy.vcf.gz",
    params:
        outdir="results/lumpy/{projectid}/{sampleid}",
    benchmark:
        "results/performance_benchmarks/lumpy_run/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/smoove.yaml"
    threads: 1
    resources:
        mem_mb="8000",
        qname="large",
    shell:
        "smoove call --outdir {params.outdir} --exclude {input.bed} --name {wildcards.sampleid} --fasta {input.fasta} -p 1 --genotype {input.bam} && "
        "mv {params.outdir}/{wildcards.sampleid}-smoove.genotyped.vcf.gz {output.vcf}"
