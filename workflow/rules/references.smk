rule download_reference_data:
    """
    Conditionally acquire reference data files from either remote
    source (e.g. S3) or somewhere else on local filesystem.

    All the interesting stuff happens in the input mapping function,
    which uses the path contents between 'reference_data' and the actual
    filename to determine which configured reference file to pull.
    """
    output:
        "reference_data/{reference_file}",
    params:
        lambda wildcards: tc.map_reference_file(wildcards, config),
    conda:
        "../envs/awscli.yaml"
    threads: 1
    resources:
        h_vmem="2000",
        qname="small",
        tmpdir="temp/",
    shell:
        'if [[ "{params}" == s3://* ]] ; then aws s3 cp {input} {output} ; else cp {input} {output} ; fi'


rule index_vcf:
    """
    Use tabix to generate tbi file for vcf.gz input
    """
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        h_vmem="2000",
        qname="small",
    shell:
        "tabix -p vcf {input}"
