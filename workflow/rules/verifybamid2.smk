rule estimate_contamination:
    """
    Use verifybamid2 to estimate contamination in samples
    """
    input:
        bam="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        fasta=config["references"]["grch37"]["fasta"],
        db_files=expand(
            "{prefix}.{suffix}",
            prefix=config["verifybamid2"]["grch37"]["db-prefix"],
            suffix=["V", "UD", "mu", "bed"],
        ),
    output:
        selfSM="{pathprefix}/contamination/{fileprefix}.vb2.selfSM",
        tsv="{pathprefix}/contamination/{fileprefix}.vb2.tsv",
    params:
        outprefix="{pathprefix}/contamination/{fileprefix}.vb2",
        db_prefix=config["verifybamid2"]["grch37"]["db-prefix"],
    conda:
        "../envs/verifybamid2.yaml"
    threads: 4
    resources:
        h_vmem="8000",
        qname="small",
    shell:
        "verifybamid2 "
        "--BamFile {input.bam} "
        "--Output {params.outprefix} "
        "--Reference {input.fasta} "
        "--SVDPrefix {params.db_prefix} "
        "--Verbose "
        "--NumThread {threads}"
