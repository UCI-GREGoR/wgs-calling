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
        r1_fastq=lambda wildcards: "{}_1_fastp.fastq".format(
            map_fastqs_to_sample_id(
                wildcards,
                manifest,
                "results/fastqs/{}/{}_R1_{}.fastq.gz".format(
                    wildcards.projectid, wildcards.prefix, wildcards.suffix
                ),
            )
        ),
        r2_fastq=lambda wildcards: "{}_2_fastp.fastq".format(
            map_fastqs_to_sample_id(
                wildcards,
                manifest,
                "results/fastqs/{}/{}_R1_{}.fastq.gz".format(
                    wildcards.projectid, wildcards.prefix, wildcards.suffix
                ),
            )
        ),
    params:
        outdir="results/fastp/{projectid}",
        sampleid=lambda wildcards: map_fastqs_to_sample_id(
            wildcards,
            manifest,
            "results/fastqs/{}/{}_R1_{}.fastq.gz".format(
                wildcards.projectid, wildcards.prefix, wildcards.suffix
            ),
        ),
        quality=10,
    conda:
        "../envs/fastp.yaml"
    threads: 1
    resources:
        h_vmem="4000",
        mem_mb="4000",
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
        "--dont_overwrite "
        "--verbose "
        "-D --overrepresentation_analysis "
        "--overrepresentation_sampling 100 "
        "--low_complexity_filter "
        "--reads_to_process 100000000 "
        "--dup_calc_accuracy 6 "
        "--trim_tail1=1 "
        "-j {params.sample}_fastp.json -h {params.sample}_fastp.html "
        "-w {threads} -z 1"
