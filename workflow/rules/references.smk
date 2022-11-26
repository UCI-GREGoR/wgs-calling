rule download_reference_data:
    """
    Conditionally acquire reference data files from either remote
    source (e.g. S3) or somewhere else on local filesystem.

    All the interesting stuff happens in the input mapping function,
    which uses the path contents between 'reference_data' and the actual
    filename to determine which configured reference file to pull.
    """
    input:
        lambda wildcards: tc.map_reference_file(wildcards, config),
    output:
        "reference_data/{reference_file}",
    benchmark:
        "results/performance_benchmarks/download_reference_data/{reference_file}.tsv"
    conda:
        "../envs/awscli.yaml"
    threads: 1
    resources:
        h_vmem="2000",
        qname="small",
        tmpdir="temp/",
    shell:
        "cp {input} {output}"


## Note that I'm switching back and forth between simple bash handling and Snakemake remote providers,
## as I'm trying to decide which one I want. The simple bash handler is as follows:
## 'if [[ "{params}" == s3://* ]] ; then aws s3 cp {params} {output} ; else cp {params} {output} ; fi'
## Also, for the bash version, the remote file needs to be a params entry; as a remote provider, it is input


rule index_vcf:
    """
    Use tabix to generate tbi file for vcf.gz input
    """
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    benchmark:
        "results/performance_benchmarks/index_vcf/{prefix}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        h_vmem="2000",
        qname="small",
    shell:
        "tabix -p vcf {input}"
