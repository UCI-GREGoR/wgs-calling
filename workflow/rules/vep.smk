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
        gz=temp("results/{prefix}.vcf.vep-annotated.gz"),
        html=temp("results/{prefix}.vcf.vep-annotated.gz_summary.html"),
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
        "{params.vep_prefix}/ensembl-vep/vep --input_file {input.vcf} --output_file {output.gz} --force_overwrite "
        "--compress_output=gzip --cache --dir_cache resources/vep/{params.reference_build} --offline "
        "--assembly={params.reference_build} --check_existing"


rule vep_format_annotation_file:
    """
    Report annotations from VEP in bcftools-friendly format
    """
    input:
        "results/{prefix}.vcf.vep-annotated.gz",
    output:
        tsv=temp("results/{prefix}.bcftools-annotation-source.tsv.gz"),
        tbi=temp("results/{prefix}.bcftools-annotation-source.tsv.gz.tbi"),
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "gunzip -c {input} | cut -f 1,13 | sed 's/_/\\t/g ; s/\\//\\t/g ; s/,/\\t/g' | "
        "awk '{{OFS = \"\\t\" ; for (i = 5 ; i <= NF ; i++) print $1,$2,$3,$4,$i}}' | "
        "awk '/\\trs/' | bgzip -c > {output.tsv} && tabix -s 1 -b 2 -e 2 {output.tsv}"


rule annotate_rsids:
    """
    Modify vep vcf outputs to place variant IDs in marker ID field
    """
    input:
        vcf="results/{prefix}.vcf.gz",
        annotations="results/{prefix}.bcftools-annotation-source.tsv.gz",
        tbi="results/{prefix}.bcftools-annotation-source.tsv.gz.tbi",
    output:
        "results/{prefix}.annotated.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "bcftools annotate -a {input.annotations} -c CHROM,POS,REF,ALT,ID -O z -o {output} {input.vcf}"
