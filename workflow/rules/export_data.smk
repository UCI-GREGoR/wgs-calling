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


def construct_export_files(suffix: str) -> list:
    """
    Use checkpoint output of linker generation rule to
    determine what files need to be constructed for a
    data export
    """
    res = []
    linker_fn = checkpoints.generate_linker.get().output[0]
    linker_df = pd.read_csv(linker_fn, sep="\t")
    res = expand(
        "results/export/{file_prefix}.{file_suffix}",
        file_prefix=linker_df["output"],
        file_suffix=suffix,
    )
    return res


rule create_bam_export:
    """
    Take bqsr bamfile and turn it into something to release
    """
    input:
        "results/bqsr/{projectid}/{sampleid}.bam",
    output:
        "results/export/{projectid}/{sampleid}_{lsid}_{sqid}.bam",
    params:
        pipeline_version=pipeline_version,
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_bam_export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "samtools reheader -c 'sed \"s/SM:{wildcards.sampleid}/SM:{params.exportid}/ ; "
        "s/LB:{wildcards.sampleid}/LB:{params.exportid}/ ; "
        "s/PU:{wildcards.sampleid}/PU:{params.exportid}/ ; "
        "\\$a@CO wgs-pipelineVersion={params.pipeline_version}\"' {input} > {output}"


rule create_snv_vcf_export:
    """
    Take snv vcf output and turn it into something to release
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sampleid}}.sorted.vcf.gz",
            toolname=config["behaviors"]["snv-caller"],
        ),
    output:
        "results/export/{projectid}/{sampleid}_{lsid}_{sqid}.vcf.gz",
    params:
        pipeline_version=pipeline_version,
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_snv_vcf_export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        'bcftools annotate -h <(echo "## wgs-pipelineVersion={params.pipeline_version}") -O v {input} | '
        'bcftools reheader -s <(echo "{wildcards.sampleid}\\t{params.exportid}") -o {output}'


rule create_export_manifest:
    """
    For export tracking and checkpointing purposes,
    report the expected fileset for a data release
    """
    input:
        bam=lambda wildcards: construct_export_files("bam"),
        bai=lambda wildcards: construct_export_files("bam.bai"),
        vcf=lambda wildcards: construct_export_files("snv.vcf.gz"),
        tbi=lambda wildcards: construct_export_files("snv.vcf.gz.tbi"),
    output:
        "results/export/manifest.tsv",
    shell:
        "echo {input.bam} {input.bai} {input.vcf} {input.tbi} > {output}"
