localrules:
    multiqc_link_ids,
    multiqc_link_ids_index_sortorder,


rule multiqc_link_ids:
    """
    Link SQ IDs to subject IDs for user convenience in multiqc reports
    """
    input:
        "results/export/linker.tsv",
    output:
        "results/multiqc/{projectid}/linker.tsv",
    shell:
        'awk -F"\\t" \'/{wildcards.projectid}/ {{print $4"\\t"$1}}\' {input} > {output}'


rule multiqc_link_ids_index_sortorder:
    """
    Link SQ IDs to subject IDs for user convenience in multiqc reports,
    but due to user feedback, have it be sorted by SQ index code instead
    of study ID.
    """
    input:
        "results/export/linker.tsv",
    output:
        "results/multiqc/{projectid}/linker_index_sortorder.tsv",
    shell:
        'awk -F"\\t" \'/{wildcards.projectid}/ {{print $4"\\t"$4"_"$1}}\' {input} > {output}'


rule run_multiqc_fastq:
    """
    Run multiqc on fastqc and fastp output for input fastqs
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc", "001_fastqc", True
        ),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc_posttrimming", "fastp_fastqc", True
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
        multiqc_config=config["multiqc-read-config"],
        id_linker="results/multiqc/{projectid}/linker.tsv",
    output:
        html="results/multiqc/{projectid}/multiqc.lane-specific.{projectid}.fastq.html",
        data_zip="results/multiqc/{projectid}/multiqc.lane-specific.{projectid}.fastq_data.zip",
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
        "../envs/multiqc.yaml" if not use_containers else None
    container:
        "docker://ewels/multiqc:v1.14" if use_containers else None
    threads: config_resources["multiqc"]["threads"]
    resources:
        mem_mb=config_resources["multiqc"]["memory"],
        qname=rc.select_queue(config_resources["multiqc"]["queue"], config_resources["queues"]),
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "--replace-names {input.id_linker} "
        "-m fastqc -m fastp "
        "-x '*.fastq.gz' -x '*.fastq' "
        "--profile-runtime --zip-data-dir "
        "-f -i 'MultiQC for Pre-alignment Data' "
        "-n {output.html}"


use rule run_multiqc_fastq as run_multiqc_fastq_index_sortorder with:
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc", "001_fastqc", True
        ),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc_posttrimming", "fastp_fastqc", True
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
        multiqc_config=config["multiqc-read-config"],
        id_linker="results/multiqc/{projectid}/linker_index_sortorder.tsv",
    output:
        html="results/multiqc/{projectid}/multiqc.lane-specific.{projectid}.fastq.html",
        data_zip="results/multiqc/{projectid}/multiqc.lane-specific.{projectid}.fastq_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_fastq_index_sortorder/{projectid}.tsv"


rule run_multiqc_alignment:
    """
    Run multiqc on all steps up to but not including variant calling
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc_combined", "fastqc", False
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc_posttrimming_combined", "fastqc", False
        ),
        verify=lambda wildcards: tc.construct_contamination_targets(wildcards, manifest),
        alignstats=tc.construct_combined_alignstats_targets,
        somalier=tc.construct_somalier_relate_targets,
        picard=lambda wildcards: tc.construct_picard_qc_targets(wildcards, manifest),
        mosdepth=lambda wildcards: tc.construct_mosdepth_targets(wildcards, manifest),
        multiqc_config=config["multiqc-alignment-config"],
        id_linker="results/multiqc/{projectid}/linker.tsv",
    output:
        html="results/multiqc/{projectid}/multiqc.combined-lanes.{projectid}.alignment.html",
        data_zip="results/multiqc/{projectid}/multiqc.combined-lanes.{projectid}.alignment_data.zip",
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
        "../envs/multiqc.yaml" if not use_containers else None
    container:
        "docker://ewels/multiqc:v1.14" if use_containers else None
    threads: config_resources["multiqc"]["threads"]
    resources:
        mem_mb=config_resources["multiqc"]["memory"],
        qname=rc.select_queue(config_resources["multiqc"]["queue"], config_resources["queues"]),
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "--replace-names {input.id_linker} "
        "-m fastqc -m fastp -m verifybamid -m picard -m somalier -m mosdepth -m custom_content "
        "--interactive "
        "-x '*.js' "
        "-x '*.bam' "
        "-x '*.fastq.gz' "
        "-x '*.pyc' "
        "-f -i 'MultiQC for Read Alignment' "
        "-n {output.html} "
        "--profile-runtime --zip-data-dir"


use rule run_multiqc_alignment as run_multiqc_alignment_index_sortorder with:
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc_combined", "fastqc", False
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, "results/fastqc_posttrimming_combined", "fastqc", False
        ),
        verify=lambda wildcards: tc.construct_contamination_targets(wildcards, manifest),
        alignstats=tc.construct_combined_alignstats_targets,
        somalier=tc.construct_somalier_relate_targets,
        picard=lambda wildcards: tc.construct_picard_qc_targets(wildcards, manifest),
        mosdepth=lambda wildcards: tc.construct_mosdepth_targets(wildcards, manifest),
        multiqc_config=config["multiqc-alignment-config"],
        id_linker="results/multiqc/{projectid}/linker_index_sortorder.tsv",
    output:
        html="results/multiqc/{projectid}/multiqc.combined-lanes.{projectid}.alignment.html",
        data_zip="results/multiqc/{projectid}/multiqc.combined-lanes.{projectid}.alignment_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_alignment_index_sortorder/{projectid}.tsv"
