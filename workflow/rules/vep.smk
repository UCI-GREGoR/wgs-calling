rule vep_download_databases:
    """
    Install backend annotation files for vep

    Note that there are some issues here since ${CONDA_PREFIX} is not available during dag build
    """
    output:
        "results/vep/{}/.download.tracker.txt".format(reference_build),
    params:
        reference_build=sm.format_reference_build(reference_build),
        vep_prefix="/opt/vep/src",
    benchmark:
        "results/performance_benchmarks/vep_download_databases/out.tsv"
    container:
        "docker://ensemblorg/ensembl-vep"
    threads: 1
    resources:
        mem_mb="1000",
        qname="small",
    shell:
        "mkdir -p resources/vep && "
        "perl {params.vep_prefix}/ensembl-vep/INSTALL.pl -a cfp -s homo_sapiens -y {params.reference_build} -c resources/vep/{params.reference_build} -g all && "
        "touch {output}"


rule vep_convert_cache:
    """
    Use tabix and a mysterious perl script to accelerate annotation
    """
    input:
        "results/vep/{}/.download.tracker.txt".format(reference_build),
    output:
        "results/vep/{}/.convert.tracker.txt".format(reference_build),
    params:
        reference_build=sm.format_reference_build(reference_build),
        vep_prefix="/opt/vep/src",
    benchmark:
        "results/performance_benchmarks/vep_convert_cache/out.tsv"
    container:
        "docker://ensemblorg/ensembl-vep"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "perl {params.vep_prefix}/ensembl-vep/convert_cache.pl -species homo_sapiens --dir resources/vep/{params.reference_build} -version all && touch {output}"


rule vep_annotate:
    """
    Run vep on a vcf file
    """
    input:
        tracker="results/vep/{}/.convert.tracker.txt".format(reference_build),
        vcf="results/{prefix}.vcf.gz",
        tbi="results/{prefix}.vcf.gz.tbi",
    output:
        temp("results/{prefix}.vcf.vep-annotated.gz"),
    params:
        reference_build=sm.format_reference_build(reference_build),
        vep_prefix="/opt/vep/src",
    container:
        "docker://ensemblorg/ensembl-vep"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "{params.vep_prefix}/ensembl-vep/vep --input_file {input.vcf} --output_file {output} --force_overwrite "
        "--compress_output=gzip --cache --dir_cache resources/vep/{params.reference_build} --offline "
        "--assembly={params.reference_build} --check_existing"


rule annotate_rsids:
    """
    Modify vep vcf outputs to place variant IDs in marker ID field
    """
    input:
        "results/{prefix}.vcf.vep-annotated.gz",
    output:
        "results/{prefix}.annotated.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "bcftools annotate ..."
