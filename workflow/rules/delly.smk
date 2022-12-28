rule delly_run:
    """
    Run delly on a single sample bam
    """
    input:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        vcf="results/delly/{projectid}/{sampleid}.delly.vcf.gz",
    conda:
        "../envs/delly.yaml"
    benchmark:
        "results/performance_benchmarks/delly_run/{projectid}/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb="16000",
        qname="small",
    shell:
        "delly call -g {input.fasta} {input.bam} | bgzip -c > {output.vcf}"
