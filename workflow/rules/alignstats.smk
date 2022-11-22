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
