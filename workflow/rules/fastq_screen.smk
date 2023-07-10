rule fastq_screen_get_references:
    """
    Use --get_genomes utility function to grab a single copy of processed reference genomes
    """
    output:
        "reference_data/FastQ_Screen_Genomes/fastq_screen.conf",
    params:
        outdir="reference_data",
    benchmark:
        "results/performance_benchmarks/fastq_screen_get_references/metrics.tsv"
    conda:
        "../envs/fastq_screen.yaml" if not use_containers else None
    container:
        "{}/fastq_screen.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["fastq_screen"]["threads"]
    resources:
        mem_mb=config_resources["fastq_screen"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["fastq_screen"]["queue"], config_resources["queues"]
        ),
    shell:
        "fastq_screen --threads {threads} --get_genomes --outdir {params.outdir} && "
        "sed -r 's|(DATABASE\\t[^\\t]+\\t).*(reference_data/.*)|\\1\\2|' {params.outdir}/fastq_screen.conf > {params.outdir}/fastq_screen.conf.tmp && "
        "mv {params.outdir}/fastq_screen.conf.tmp {params.outdir}/fastq_screen.conf"


rule fastq_screen_run:
    """
    Run fastq-screen on a single fastq to estimate species contributions
    """
    input:
        fastq="results/fastqs/{projectid}/{sampleid}_{lane}_{read}_001.fastq.gz",
        config="reference_data/FastQ_Screen_Genomes/fastq_screen.conf",
    output:
        expand(
            "results/fastq_screen/{{projectid}}/{{sampleid}}_{{lane}}_{{read}}_001_screen.{suffix}",
            suffix=["txt", "png", "html"],
        ),
    params:
        outdir="results/fastq_screen/{projectid}",
    conda:
        "../envs/fastq_screen.yaml" if not use_containers else None
    container:
        "{}/fastq_screen.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["fastq_screen"]["threads"]
    resources:
        mem_mb=config_resources["fastq_screen"]["memory"],
        qname=rc.select_queue(config_resources["fastq_screen"]["queue"], config_resources["queues"]),
    shell:
        "fastq_screen --threads {threads} --conf {input.config} --aligner bowtie2 --outdir {params.outdir} {input.fastq}"