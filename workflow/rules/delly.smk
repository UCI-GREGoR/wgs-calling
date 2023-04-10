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
    container:
        "docker://dellytools/delly:latest"
    benchmark:
        "results/performance_benchmarks/delly_run/{projectid}/{sampleid}.tsv"
    threads: config_resources["delly"]["threads"]
    resources:
        mem_mb=config_resources["delly"]["memory"],
        qname=rc.select_queue(config_resources["delly"]["queue"], config_resources["queues"]),
    shell:
        "delly call -g {input.fasta} -x {input.bed} {input.bam} | bcftools view -i 'F_MISSING<0.5' -O z -o {output.vcf}"
