rule run_alignstats:
    """
    Run alignstats utility on post-markdups bams
    """
    input:
        bam="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam.bai",
    output:
        json="{pathprefix}/alignstats/{fileprefix}.bwa2a.alignstats.json",
    params:
        min_qual=20,
    conda:
        "../envs/alignstats.yaml"
    threads: 4
    resources:
        h_vmem="8000",
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
    threads: 1
    resources:
        h_vmem="2000",
        qname="small",
    script:
        "../scripts/alignstats_json_to_tsv.py"
