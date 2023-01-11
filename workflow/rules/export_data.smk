checkpoint generate_linker:
    """
    From a sample logbook, generate a simple linker
    between various sample ID types
    """
    input:
        logbook=config["sample-logbook"],
    output:
        linker="results/export/linker.tsv",
    benchmark:
        "results/performance_benchmarks/generate_linker/linker.tsv"
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    script:
        "../scripts/construct_linker_from_labbook.R"
