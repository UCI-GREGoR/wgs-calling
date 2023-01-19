rule run_multiqc_fastq:
    """
    Run multiqc on fastqc and fastp output for input fastqs
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(wildcards, manifest),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_posttrimming_targets(
            wildcards, manifest
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
        multiqc_config=config["multiqc-read-config"],
    output:
        html="results/multiqc/{projectid}/multiqc.fastq.html",
        data_zip="results/multiqc/{projectid}/multiqc.fastq_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_fastq/{projectid}.tsv"
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{{projectid}}",
                    toolname=["fastqc", "fastp", "fastqc_posttrimming"],
                )
            )
        ),
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "-m fastqc -m fastp "
        "-x '*.fastq.gz' -x '*.fastq' "
        "--profile-runtime --zip-data-dir "
        "-f -i 'MultiQC for Pre-alignment Data' "
        "-n {output.html}"


rule run_multiqc_alignment:
    """
    Run multiqc on all steps up to but not including variant calling
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_combined_targets(wildcards, manifest),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_posttrimming_combined_targets(
            wildcards, manifest
        ),
        verify=lambda wildcards: tc.construct_contamination_targets(wildcards, manifest),
        alignstats=tc.construct_combined_alignstats_targets,
        somalier=tc.construct_somalier_relate_targets,
        picard=lambda wildcards: tc.construct_picard_qc_targets(wildcards, manifest),
        mosdepth=lambda wildcards: tc.construct_mosdepth_targets(wildcards, manifest),
        multiqc_config=config["multiqc-alignment-config"],
    output:
        html="results/multiqc/{projectid}/multiqc.alignment.html",
        data_zip="results/multiqc/{projectid}/multiqc.alignment_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_alignment/{projectid}.tsv"
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{{projectid}}",
                    toolname=[
                        "fastqc_combined",
                        "fastqc_posttrimming_combined",
                        "collectmultiplemetrics",
                        "collectgcbiasmetrics",
                        "collectwgsmetrics",
                        "somalier",
                        "mosdepth",
                        "contamination",
                        "alignstats",
                        "markdups",
                    ],
                )
            )
        ),
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb="4000",
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


rule run_multiqc_calling:
    """
    Run multiqc on variant calling only
    """
    input:
        vcf_raw=tc.construct_snv_targets(config, manifest),
        vcf_export=lambda wildcards: ed.construct_export_files(wildcards, manifest, "snv.vcf.stats"),
        vcf_nonexport=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, "snv.vcf.stats"
        ),
        multiqc_config=config["multiqc-calling-config"],
    output:
        html="results/multiqc/{projectid}/multiqc.calling.html",
        data_zip="results/multiqc/{projectid}/multiqc.calling_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_calling/{projectid}.tsv"
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{{projectid}}",
                    toolname=[
                        "export",
                        "nonexport",
                    ],
                )
            )
        ),
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "-m bcftools/stats "
        "--interactive "
        "-f -i 'MultiQC for Variant Calling' "
        "-n {output.html} "
        "--profile-runtime --zip-data-dir"
