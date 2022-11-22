rule estimate_contamination:
    """
    Use verifybamid2 to estimate contamination in samples
    """
    input:
        bam="{pathprefix}/contamination/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/contamination/{fileprefix}.mrkdup.sort.bam.bai",
        fasta=config["references"]["grch37"]["fasta"],
        db_prefix=config["verifybamid"]["grch37"]["db-prefix"],
    output:
        selfSM="{pathprefix}/contamination/{fileprefix}.vb2.selfSM",
        tsv="{pathprefix}/contamination/{fileprefix}.vb2.tsv",
    params:
        outprefix="{pathprefix}/contamination/{fileprefix}.vb2",
    conda:
        "../envs/verifybamid2.yaml"
    threads: 4
    resources:
        h_vmem="8000",
        qname="small",
    shell:
        "VerifyBamID "
        "--BamFile {input.bam} "
        "--Output {params.outprefix} "
        "--Reference {input.fasta} "
        "--SVDPrefix {input.db_prefix} "
        "--Verbose "
        "--NumThread {threads}"
