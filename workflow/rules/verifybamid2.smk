rule estimate_contamination:
    """
    Use verifybamid2 to estimate contamination in samples
    """
    input:
        bam="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{fileprefix}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/fasta".format(reference_build),
        db_files=expand(
            "reference_data/verifybamid2/{genome}/db.{suffix}",
            genome=reference_build,
            suffix=["V", "UD", "mu", "bed"],
        ),
    output:
        selfSM="{pathprefix}/contamination/{fileprefix}.vb2.selfSM",
        ancestry="{pathprefix}/contamination/{fileprefix}.vb2.Ancestry",
    params:
        outprefix="{pathprefix}/contamination/{fileprefix}.vb2",
        db_prefix="reference_data/verifybamid2/{}/db",
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
        "--DisableSanityCheck "
        "--NumThread {threads}"
