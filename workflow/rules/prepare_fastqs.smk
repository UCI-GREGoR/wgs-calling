rule sort_input_bam:
    """
    For input files provided as pre-aligned bams: the expectation is that these
    bams will be aligned to the wrong genome, and need to be converted into fastqs
    in preparation for re-alignment. As a preliminary requirement, sort the bam.

    This bam will be passed along to samtools fixmate, and as such needs to be
    sorted with -n.
    """
    input:
        bam=lambda wildcards: tc.locate_input_bam(wildcards, manifest, True),
    output:
        bam=temp("results/input_bams/{projectid}/{sampleid}.sorted.bam"),
    benchmark:
        "results/performance_benchmarks/sort_input_bam/{projectid}/{sampleid}.tsv"
    params:
        tmpdir=tempDir,
        sort_m="{}M".format(
            int(
                float(config_resources["samtools"]["memory"])
                / (2 * float(config_resources["samtools"]["threads"]))
            )
        ),
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/bwa.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools_sort"]["threads"]
    resources:
        mem_mb=config_resources["samtools_sort"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools_sort"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        "samtools sort -@ {threads} -T {params.tmpdir} -m {params.sort_m} -n -o {output.bam} -O bam {input.bam}"


rule fix_mate_bam:
    """
    Some, but not all, of the bams this pipeline has been receiving require samtools
    fixmate to be run on them before running samtools fastq; without this step,
    the vast majority of reads are removed as singletons, in spite of the fact
    that samtools flagstat doesn't think that singleton behavior should happen.
    """
    input:
        bam="results/input_bams/{projectid}/{sampleid}.sorted.bam",
    output:
        bam=temp("results/input_bams/{projectid}/{sampleid}.fixmate.bam"),
    benchmark:
        "results/performance_benchmarks/fix_mate_bam/{projectid}/{sampleid}.tsv"
    params:
        tmpdir=tempDir,
        sort_m="{}M".format(
            int(
                float(config_resources["samtools"]["memory"])
                / (2 * float(config_resources["samtools"]["threads"]))
            )
        ),
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/bwa.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools"]["threads"]
    resources:
        mem_mb=config_resources["samtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools"]["queue"], config_resources["queues"]
        ),
        tmpdir=tempDir,
    shell:
        "samtools fixmate -@ 16 {input.bam} {output.bam}"


checkpoint input_bam_sample_lanes:
    """
    For input files provided as pre-aligned bams: the expectation is that these
    bams will be aligned to the wrong genome, and need to be converted into fastqs
    in preparation for re-alignment. To determine the expected files after splitting
    by lane, sniff the bam for read names and determine which lanes are reportedly present.
    """
    input:
        "results/input_bams/{projectid}/{sampleid}.fixmate.bam",
    output:
        "results/fastqs_from_bam/{projectid}/{sampleid}_expected-lanes.tsv",
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
        'samtools view {input} | cut -f 4 -d ":" | sort | uniq > {output}'


rule input_bam_to_split_fastq:
    """
    For input files provided as pre-aligned bams: the expectation is that these
    bams will be aligned to the wrong genome, and need to be converted into fastqs
    in preparation for re-alignment. After sorting the bam, convert it to fastq,
    split by lane, bgzip compressed.
    """
    input:
        bam="results/input_bams/{projectid}/{sampleid}.fixmate.bam",
        expected="results/fastqs_from_bam/{projectid}/{sampleid}_expected-lanes.tsv",
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
        "samtools fastq -@ {threads} -s /dev/null -{params.off_target_read_flag} /dev/null -0 /dev/null -n {input.bam} | "
        'awk -v target={wildcards.lane} \'BEGIN {{FS = ":"}} {{lane = $4 ; if (lane == target) {{print}} ; '
        "for (i = 1 ; i <= 3 ; i++) {{getline ; if (lane == target) {{print}}}}}}' | "
        "bgzip -c > {output}"


use rule copy_fastqs as copy_combined_fastqs with:
    output:
        fastq=temp(
            "results/imported_fastqs/{projectid}/{sampleid}_{lane}_R{readgroup}_{suffix}.fastq.gz"
        ),
    benchmark:
        "results/performance_benchmarks/copy_combined_fastqs/{projectid}/{sampleid}_{lane}_R{readgroup}_{suffix}.fastq.tsv"


checkpoint input_fastq_sample_lanes:
    """
    For input files provided as combined fastqs: the expectation is that these
    fastqs will need to be split into per-lane fastqs for both improved performance
    and finer-grained QC. To determine the expected files after splitting
    by lane, sniff the fastq for read names and determine which lanes are reportedly present.
    """
    input:
        "results/imported_fastqs/{projectid}/{sampleid}_combined_{readgroup}_001.fastq.gz",
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
        "results/imported_fastqs/{projectid}/{sampleid}_combined_{readgroup}_001.fastq.gz",
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
