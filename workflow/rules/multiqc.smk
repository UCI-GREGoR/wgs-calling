localrules:
    multiqc_link_ids_index_sortorder,


rule multiqc_link_ids_index_sortorder:
    """
    Link original sequencing index IDs to subject IDs for user convenience in multiqc reports,
    but due to user feedback, have it be sorted by sequencing index code instead
    of study ID.
    """
    input:
        "results/export/linker.tsv",
    output:
        "results/multiqc/{projectid}/linker_index_sortorder.tsv",
    shell:
        'awk -F"\\t" \'/{wildcards.projectid}/ && $1 != $4 {{print $4"\\t"$4"_"$1}}\' {input} > {output}'


rule run_multiqc_fastq_lane_specific:
    """
    Run multiqc on fastqc and fastp output for input fastqs split by lane
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, checkpoints, "results/fastqc", "001_fastqc", True
        ),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, checkpoints, "results/fastqc_posttrimming", "fastp_fastqc", True
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest, checkpoints),
        fastq_screen=lambda wildcards: tc.construct_fastq_screen_targets(
            wildcards, manifest, checkpoints, False
        ),
        multiqc_config="config/multiqc_read_config_lane_specific.yaml",
        id_linker="results/multiqc/{projectid}/linker_index_sortorder.tsv",
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
                    toolname=["fastqc", "fastp", "fastqc_posttrimming", "fastq_screen"],
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
        qname=lambda wildcards: rc.select_queue(
            config_resources["multiqc"]["queue"], config_resources["queues"]
        ),
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "--replace-names {input.id_linker} "
        "-m fastqc -m fastp -m fastq_screen "
        "-x '*.fastq.gz' -x '*.fastq' "
        "--profile-runtime --zip-data-dir "
        "-f -i 'MultiQC for Pre-alignment Data' "
        "-n {output.html}"


use rule run_multiqc_fastq_lane_specific as run_multiqc_fastq_combined_lanes with:
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, checkpoints, "results/fastqc_combined", "001_fastqc", False
        ),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards,
            manifest,
            checkpoints,
            "results/fastqc_posttrimming_combined",
            "001_fastqc",
            False,
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest, checkpoints),
        fastq_screen=lambda wildcards: tc.construct_fastq_screen_targets(
            wildcards, manifest, checkpoints, True
        ),
        multiqc_config="config/multiqc_read_config_combined_lanes.yaml",
        id_linker="results/multiqc/{projectid}/linker_index_sortorder.tsv",
    output:
        html="results/multiqc/{projectid}/multiqc.combined-lanes.{projectid}.fastq.html",
        data_zip="results/multiqc/{projectid}/multiqc.combined-lanes.{projectid}.fastq_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_fastq_combined_lanes/{projectid}.tsv"
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{{projectid}}",
                    toolname=[
                        "fastqc_combined",
                        "fastp",
                        "fastqc_posttrimming_combined",
                        "fastq_screen_combined",
                    ],
                )
            )
        ),


rule run_multiqc_alignment_combined_lanes:
    """
    Run multiqc on all steps up to but not including variant calling, with all lanes combined
    """
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, checkpoints, "results/fastqc_combined", "001_fastqc", False
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest, checkpoints),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards,
            manifest,
            checkpoints,
            "results/fastqc_posttrimming_combined",
            "001_fastqc",
            False,
        ),
        verify=lambda wildcards: tc.construct_contamination_targets(wildcards, manifest),
        alignstats=tc.construct_combined_alignstats_targets,
        somalier=tc.construct_somalier_relate_targets,
        picard=lambda wildcards: tc.construct_picard_qc_targets(wildcards, manifest),
        mosdepth=lambda wildcards: tc.construct_mosdepth_targets(wildcards, manifest),
        fastq_screen=lambda wildcards: tc.construct_fastq_screen_targets(
            wildcards, manifest, checkpoints, True
        ),
        multiqc_config="config/multiqc_alignment_config_combined_lanes.yaml",
        id_linker="results/multiqc/{projectid}/linker_index_sortorder.tsv",
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
                        "fastq_screen_combined",
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
        qname=lambda wildcards: rc.select_queue(
            config_resources["multiqc"]["queue"], config_resources["queues"]
        ),
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "--replace-names {input.id_linker} "
        "-m fastqc -m fastp -m verifybamid -m picard -m somalier -m mosdepth -m fastq_screen -m custom_content "
        "--interactive "
        "-x '*.js' "
        "-x '*.bam' "
        "-x '*.fastq.gz' "
        "-x '*.pyc' "
        "-f -i 'MultiQC for Read Alignment' "
        "-n {output.html} "
        "--profile-runtime --zip-data-dir"


use rule run_multiqc_alignment_combined_lanes as run_multiqc_alignment_lane_specific with:
    input:
        fastqc=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, checkpoints, "results/fastqc", "001_fastqc", True
        ),
        fastp=lambda wildcards: tc.construct_fastp_targets(wildcards, manifest, checkpoints),
        fastqc_posttrimming=lambda wildcards: tc.construct_fastqc_targets(
            wildcards, manifest, checkpoints, "results/fastqc_posttrimming", "fastp_fastqc", True
        ),
        verify=lambda wildcards: tc.construct_contamination_targets(wildcards, manifest),
        alignstats=tc.construct_combined_alignstats_targets,
        somalier=tc.construct_somalier_relate_targets,
        picard=lambda wildcards: tc.construct_picard_qc_targets(wildcards, manifest),
        mosdepth=lambda wildcards: tc.construct_mosdepth_targets(wildcards, manifest),
        fastq_screen=lambda wildcards: tc.construct_fastq_screen_targets(
            wildcards, manifest, checkpoints, False
        ),
        multiqc_config="config/multiqc_alignment_config_lane_specific.yaml",
        id_linker="results/multiqc/{projectid}/linker_index_sortorder.tsv",
    output:
        html="results/multiqc/{projectid}/multiqc.lane-specific.{projectid}.alignment.html",
        data_zip="results/multiqc/{projectid}/multiqc.lane-specific.{projectid}.alignment_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_alignment_lane_specific/{projectid}.tsv"
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{{projectid}}",
                    toolname=[
                        "fastqc",
                        "fastqc_posttrimming",
                        "collectmultiplemetrics",
                        "collectgcbiasmetrics",
                        "collectwgsmetrics",
                        "somalier",
                        "mosdepth",
                        "fastq_screen",
                        "contamination",
                        "alignstats",
                        "markdups",
                    ],
                )
            )
        ),
