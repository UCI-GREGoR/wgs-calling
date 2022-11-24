rule run_multiqc_fastq:
    """
    Run multiqc on fastqc and fastp output for input fastqs
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(wildcards, manifest),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
    output:
        html="results/multiqc/{projectid}/multiqc.fastq.html",
        data_zip="results/multiqc/{projectid}/multiqc.fastq_data.zip",
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{{projectid}}",
                    toolname=["fastqc", "fastp"],
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
        "multiqc {params.target_dirs} "
        "-m fastqc -m fastp "
        "-x '*.fastq.gz' -x '*.fastq' "
        "--profile-runtime --zip-data-dir "
        "-f -i 'MultiQC for Raw Fastqs' "
        "-n {output.html}"


rule run_multiqc_alignment:
    """
    Run multiqc on all steps up to but not including variant calling
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(wildcards, manifest),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
        verify=lambda wildcards: tc.construct_contamination_targets(wildcards, manifest),
        alignstats=tc.construct_combined_alignstats_targets,
        somalier=tc.construct_somalier_relate_targets,
        picard=lambda wildcards: tc.construct_picard_qc_targets(wildcards, manifest),
        mosdepth=lambda wildcards: tc.construct_mosdepth_targets(wildcards, manifest),
        multiqc_config=config["multiqc-config"],
    output:
        html="results/multiqc/{projectid}/multiqc.alignment.html",
        data_zip="results/multiqc/{projectid}/multiqc.alignment_data.zip",
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{{projectid}}",
                    toolname=[
                        "fastqc",
                        "fastp",
                        "collectmultiplemetrics",
                        "collectgcbiasmetrics",
                        "collectwgsmetrics",
                        "somalier",
                        "mosdepth",
                    ],
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
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "-m fastqc -m fastp -m verifybamid -m picard -m somalier -m mosdepth -m custom_content "
        "--interactive "
        "-x '*.js' "
        "-x '*.bam' "
        "-x '*.fastq.gz' "
        "-x '*.pyc' "
        "-f -i 'MultiQC for Read Alignment' "
        "-n {output.html} "
        "--profile-runtime --zip-data-dir"
