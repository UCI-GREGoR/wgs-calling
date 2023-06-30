"""
For input files provided as pre-aligned bams: the expectation is that these
bams will be aligned to the wrong genome, and need to be converted into fastqs
in preparation for re-alignment. As a preliminary requirement, sort the bam.
"""


use rule sort_bam as sort_input_bam with:
    input:
        bam=lambda wildcards: tc.locate_input_bam(wildcards, manifest, True),
    output:
        bam=temp("results/input_bams/{projectid}/{sampleid}.sorted.bam"),
    benchmark:
        "results/performance_benchmarks/sort_input_bam/{projectid}/{sampleid}.tsv"


checkpoint input_bam_sample_lanes:
    """
    For input files provided as pre-aligned bams: the expectation is that these
    bams will be aligned to the wrong genome, and need to be converted into fastqs
    in preparation for re-alignment. To determine the expected files after splitting
    by lane, sniff the bam for read names and determine which lanes are reportedly present.
    """
    input:
        "results/input_bams/{projectid}/{sampleid}.sorted.bam",
    output:
        temp("results/fastqs_from_bam/{projectid}/{sampleid}_expected-lanes.tsv"),
    benchmark:
        "results/performance_benchmarks/input_bam_sample_lanes/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/samtools.sif".format(apptainer_images) if use_containers else None
    threads: 1
    resources:
        mem_mb=2000,
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools"]["queue"], config_resources["queues"]
        ),
    shell:
        'samtools head -h 0 -n 100000 {input} | cut -f 4 -d ":" | sort | uniq > {output}'


rule input_bam_to_split_fastq:
    """
    For input files provided as pre-aligned bams: the expectation is that these
    bams will be aligned to the wrong genome, and need to be converted into fastqs
    in preparation for re-alignment. After sorting the bam, convert it to fastq,
    split by lane, bgzip compressed.
    """
    input:
        "results/input_bams/{projectid}/{sampleid}.sorted.bam",
    output:
        "results/fastqs_from_bam/{projectid}/{sampleid}_L00{lane}_{readgroup}_001.fastq.gz",
    benchmark:
        "results/performance_benchmarks/input_bam_to_split_fastq/{projectid}/{sampleid}_L00{lane}_{readgroup}.tsv"
    params:
        off_target_read_flag=lambda wildcards: 3 - int(wildcards.readgroup.strip("R")),
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/samtools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools"]["threads"]
    resources:
        mem_mb=config_resources["samtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "samtools fastq -@ {threads} -s /dev/null -{params.off_target_read_flag} /dev/null -0 /dev/null -n {input} | "
        'awk -v target={wildcards.lane} \'BEGIN {{FS = ":"}} {{lane = $4 ; if (lane == target) {{print}} ; '
        "for (i = 1 ; i <= 3 ; i++) {{getline ; if (lane == target) {{print}}}}}}' | "
        "bgzip -c > {output}"


checkpoint input_fastq_sample_lanes:
    """
    For input files provided as combined fastqs: the expectation is that these
    fastqs will need to be split into per-lane fastqs for both improved performance
    and finer-grained QC. To determine the expected files after splitting
    by lane, sniff the fastq for read names and determine which lanes are reportedly present.
    """
    input:
        lambda wildcards: manifest.query(
            'projectid == "{}" and sampleid == "{}"'.format(
                wildcards.projectid, wildcards.sampleid
            )
        )[wildcards.readgroup.lower()].to_list()[0],
    output:
        temp("results/fastqs_from_fastq/{projectid}/{sampleid}_{readgroup}_expected-lanes.tsv"),
    benchmark:
        "results/performance_benchmarks/input_fastq_sample_lanes/{projectid}/{sampleid}_{readgroup}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        qname=lambda wildcards: rc.select_queue("small", config_resources["queues"]),
    shell:
        "gunzip -c {input} | awk 'NF > 1 {{print $1}}' | cut -f 4 -d ':' | sort | uniq > {output}"


rule input_fastq_to_split_fastq:
    """
    For input files provided as pre-aligned bams: the expectation is that these
    bams will be aligned to the wrong genome, and need to be converted into fastqs
    in preparation for re-alignment. After sorting the bam, convert it to fastq,
    split by lane, bgzip compressed.

    Based on observations, it seems like these fastqs are minimally incompatible
    with fastp. I'm testing why, but bbtools can fix it; but that needs to be
    invoked separately.
    """
    input:
        lambda wildcards: manifest.query(
            'projectid == "{}" and sampleid == "{}"'.format(
                wildcards.projectid, wildcards.sampleid
            )
        )[wildcards.readgroup.lower()].to_list()[0],
    output:
        temp("results/bbtools_input/{projectid}/{sampleid}_L00{lane}_{readgroup}_001.fastq.gz"),
    benchmark:
        "results/performance_benchmarks/input_fastq_to_split_fastq/{projectid}/{sampleid}_L00{lane}_{readgroup}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: 1
    resources:
        mem_mb=1000,
        qname=lambda wildcards: rc.select_queue("small", config_resources["queues"]),
    shell:
        "gunzip -c {input} | "
        'awk \'BEGIN {{FS = ":"}} {{lane = $4 ; if ( lane == "{wildcards.lane}" ) {{ print }} ; '
        'for (i = 1 ; i <= 3 ; i++) {{getline ; if ( lane == "{wildcards.lane}" ) {{ print }}}}}}\' | '
        "bgzip -c > {output}"


rule bbtools_repair_fastqs:
    """
    Testing the use of bbtools repair.sh to fix incompatibilities between external combined fastqs
    and fastp.
    """
    input:
        R1="results/bbtools_input/{projectid}/{sampleid}_L00{lane}_R1_001.fastq.gz",
        R2="results/bbtools_input/{projectid}/{sampleid}_L00{lane}_R2_001.fastq.gz",
    output:
        R1="results/fastqs_from_fastq/{projectid}/{sampleid}_L00{lane}_R1_001.fastq.gz",
        R2="results/fastqs_from_fastq/{projectid}/{sampleid}_L00{lane}_R2_001.fastq.gz",
        singletons=temp(
            "results/fastqs_from_fastq/{projectid}/{sampleid}_L00{lane}_singletons_001.fastq.gz"
        ),
    benchmark:
        "results/performance_benchmarks/bbtools_repair_fastqs/{projectid}/{sampleid}_L00{lane}.tsv"
    conda:
        "../envs/bbtools.yaml" if not use_containers else None
    threads: config_resources["bbtools"]["threads"]
    resources:
        mem_mb=config_resources["bbtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bbtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "repair.sh in1={input.R1} in2={input.R2} out1={output.R1} out2={output.R2} outs={output.singletons} repair"
