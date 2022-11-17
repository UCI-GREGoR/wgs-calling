rule run_fastp:
    """
    Run fastp on a paired end library
    """
    input:
        r1="results/fastqs/{projectid}/{prefix}_R1_{suffix}.fastq.gz",
        r2="results/fastqs/{projectid}/{prefix}_R2_{suffix}.fastq.gz",
    output:
        r1_zip="results/fastp/{projectid}/{prefix}_R1_{suffix}_fastp.zip",
        r2_zip="results/fastp/{projectid}/{prefix}_R2_{suffix}_fastp.zip",
        r1_html="results/fastp/{projectid}/{prefix}_R1_{suffix}_fastp.html",
        r2_html="results/fastp/{projectid}/{prefix}_R2_{suffix}_fastp.html",
        r1_json="results/fastp/{projectid}/{prefix}_R1_{suffix}_fastp.json",
        r2_json="results/fastp/{projectid}/{prefix}_R2_{suffix}_fastp.json",
    params:
        outdir="results/fastp/{projectid}",
        sampleid=tc.map_fastqs_to_sampleid,
        r1_fastq=lambda wildcards: "{}_1_fastp.fastq".format(tc.map_fastqs_to_sampleid(wildcards)),
        r2_fastq=lambda wildcards: "{}_2_fastp.fastq".format(tc.map_fastqs_to_sampleid(wildcards)),
        quality=10,
        dup_calc_accuracy=3,
    conda:
        "../envs/fastp.yaml"
    threads: 2
    resources:
        h_vmem="16000",
        mem_mb="16000",
        qname="small",
    shell:
        "mkdir -p {params.outdir} && "
        "fastp -i {input.r1} -I {input.r2} "
        "-o {params.sampleid}_1_fastp.fastq -O {params.sampleid}_2_fastp.fastq "
        "--unpaired1 {params.sampleid}_1_unp_fastp.fastq --unpaired2 {params.sampleid}_2_unp_fastp.fastq "
        "--failed_out {params.sampleid}_failed.fastq "
        "--detect_adapter_for_pe "
        "-q {params.quality} -u 60 "
        "--trim_poly_g "
        "--verbose "
        "-D --overrepresentation_analysis "
        "--overrepresentation_sampling 100 "
        "--low_complexity_filter "
        "--reads_to_process 100000000 "
        "--dup_calc_accuracy {params.dup_calc_accuracy} "
        "--trim_tail1=1 "
        "-j {params.sampleid}_fastp.json -h {params.sampleid}_fastp.html "
        "-w {threads} -z 1"
