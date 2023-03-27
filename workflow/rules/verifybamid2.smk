rule estimate_contamination:
    """
    Use verifybamid2 to estimate contamination in samples
    """
    input:
        bam="results/bqsr/{fileprefix}.bam",
        bai="results/bqsr/{fileprefix}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        db_files=expand(
            "reference_data/verifybamid2/{genome}/ref.db.{suffix}",
            genome=reference_build,
            suffix=["V", "UD", "mu", "bed"],
        ),
    output:
        selfSM="results/contamination/{fileprefix}.vb2.selfSM",
        ancestry="results/contamination/{fileprefix}.vb2.Ancestry",
    benchmark:
        "results/performance_benchmarks/estimate_contamination/{fileprefix}.tsv"
    params:
        outprefix="results/contamination/{fileprefix}.vb2",
        db_prefix="reference_data/verifybamid2/{}/ref.db".format(reference_build),
    conda:
        "../envs/verifybamid2.yaml"
    threads: config_resources["verifybamid2"]["threads"]
    resources:
        mem_mb=config_resources["verifybamid2"]["memory"],
        qname=rc.select_queue(config_resources["verifybamid2"]["queue"]),
    shell:
        "verifybamid2 "
        "--BamFile {input.bam} "
        "--Output {params.outprefix} "
        "--Reference {input.fasta} "
        "--SVDPrefix {params.db_prefix} "
        "--Verbose "
        "--DisableSanityCheck "
        "--NumThread {threads}"
