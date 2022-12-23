rule run_alignstats:
    """
    Run alignstats utility on post-markdups bams
    """
    input:
        bam="results/bqsr/{fileprefix}.bam",
        bai="results/bqsr/{fileprefix}.bam.bai",
    output:
        json="results/alignstats/{fileprefix}.bwa2a.alignstats.json",
    benchmark:
        "results/performance_benchmarks/run_alignstats/{fileprefix}.tsv"
    params:
        min_qual=20,
    conda:
        "../envs/alignstats.yaml"
    threads: 4
    resources:
        mem_mb="8000",
        qname="small",
    shell:
        "alignstats -C -U "
        "-i {input.bam} "
        "-o {output.json} "
        "-j bam "
        "-v -P 10 -p 10 "
        "-q {params.min_qual}"


rule merge_alignstats:
    """
    Combine json-format alignstats output into a single big table
    """
    input:
        json=lambda wildcards: tc.construct_alignstats_targets(wildcards, manifest),
    output:
        tsv="results/alignstats/{projectid}/alignstats_summary_mqc.tsv",
    benchmark:
        "results/performance_benchmarks/merge_alignstats/{projectid}/alignstats_summary_mqc.tsv"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    script:
        "../scripts/alignstats_json_to_tsv.py"
