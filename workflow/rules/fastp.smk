rule run_fastp:
    """
    Run fastp on a paired end library
    """
    input:
        r1="results/fastqs/{projectid}/{sampleid}_{lane}_R1_001.fastq.gz",
        r2="results/fastqs/{projectid}/{sampleid}_{lane}_R2_001.fastq.gz",
    output:
        html="results/fastp/{projectid}/{sampleid}_{lane}_fastp.html",
        json="results/fastp/{projectid}/{sampleid}_{lane}_fastp.json",
        r1_fastq="results/fastp/{projectid}/{sampleid}_{lane}_R1_fastp.fastq",
        r2_fastq="results/fastp/{projectid}/{sampleid}_{lane}_R2_fastp.fastq",
        r1_unp_fastq="results/fastp/{projectid}/{sampleid}_{lane}_R1_unp_fastp.fastq",
        r2_unp_fastq="results/fastp/{projectid}/{sampleid}_{lane}_R2_unp_fastp.fastq",
        failed_fastq="results/fastp/{projectid}/{sampleid}_{lane}_failed.fastq",
    benchmark:
        "results/performance_benchmarks/run_fastp/{projectid}/{sampleid}_{lane}.tsv"
    params:
        outprefix="results/fastp/{projectid}/{sampleid}_{lane}",
        quality=10,
    conda:
        "../envs/fastp.yaml"
    threads: 2
    resources:
        mem_mb="8000",
        qname="small",
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-o {output.r1_fastq} -O {output.r2_fastq} "
        "--unpaired1 {output.r1_unp_fastq} --unpaired2 {output.r2_unp_fastq} "
        "--failed_out {output.failed_fastq} "
        "-q {params.quality} "
        "--trim_poly_g "
        "--verbose "
        "--overrepresentation_analysis "
        "--overrepresentation_sampling 100 "
        "-5 -3 "
        "-j {output.json} -h {output.html} "
        "-w {threads} -z 1"
