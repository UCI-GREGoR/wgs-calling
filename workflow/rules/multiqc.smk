rule run_fastq_multiqc:
    """
    Run multiqc on fastqc and fastp output for input fastqs
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(manifest),
        fastp=lambda wildcards: tc.construct_fastp_targets(manifest),
    output:
        "results/multiqc/multiqc.fastq.html",
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{projectid}",
                    toolname=["fastqc", "fastp"],
                    projectid=manifest["projectid"],
                )
            )
        ),
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    resources:
        h_vmem="4000",
        qname="small",
    shell:
        "multiqc {params.target_dirs} -p "
        "-m fastqc -m fastp "
        "-x '*.fastq.gz' -x '*.fastq' "
        "--profile-runtime --zip-data-dir "
        "-f -i 'MultiQC for Raw Fastqs' "
        "-n {output}"
