rule delly_run:
    """
    Run delly on a single sample bam
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
        bed="reference_data/delly/{}/ref.exclude.bed".format(reference_build),
    output:
        vcf="results/delly/{projectid}/{sampleid}.delly.vcf.gz",
    conda:
        "../envs/delly.yaml"
    benchmark:
        "results/performance_benchmarks/delly_run/{projectid}/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb="12000",
        qname="small",
    shell:
        "delly call -g {input.fasta} -x {input.bed} {input.bam} | bcftools view -i 'F_MISSING<0.5' -O z -o {output.vcf}"
