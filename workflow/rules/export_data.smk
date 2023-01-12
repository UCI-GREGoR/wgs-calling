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


def construct_export_files(wildcards, manifest: pd.DataFrame, suffix: str) -> list:
    """
    Use checkpoint output of linker generation rule to
    determine what files need to be constructed for a
    data export
    """
    res = []
    linker_fn = checkpoints.generate_linker.get().output[0]
    linker_df = pd.read_csv(linker_fn, sep="\t")
    subjectids = manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"].to_list()
    targets = linker_df.loc[
        (linker_df["ru"] == wildcards.projectid) & [x in subjectids for x in linker_df["sq"]],
        "output",
    ]
    res = expand(
        "results/export/{projectid}/{file_prefix}.{file_suffix}",
        projectid=wildcards.projectid,
        file_prefix=targets.to_list(),
        file_suffix=suffix,
    )
    return res


rule create_bam_export:
    """
    Take bqsr bamfile and turn it into something to release
    """
    input:
        "results/bqsr/{projectid}/{sqid}.bam",
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
        "samtools reheader -c 'sed \"s/SM:{wildcards.sqid}/SM:{params.exportid}/ ; "
        "s/LB:{wildcards.sqid}/LB:{params.exportid}/ ; "
        "s/PU:{wildcards.sqid}/PU:{params.exportid}/ ; "
        "\\$a@CO\\twgs-pipelineVersion={params.pipeline_version}\"' {input} > {output}"


rule create_bai_export:
    """
    Index export bam file

    This isn't using rule inheritance from the rule in the picard code
    due to idiosyncratic requirements of loading order of rules with inheritance
    """
    input:
        bam="results/export/{prefix}.bam",
    output:
        bai="results/export/{prefix}.bai",
    benchmark:
        "results/performance_benchmarks/create_bai_export/{prefix}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    resources:
        mem_mb="8000",
        qname="small",
    shell:
        "samtools index -@ {threads} -b -o {output.bai} {input.bam}"


rule create_snv_vcf_export:
    """
    Take snv vcf output and turn it into something to release
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sqid}}.sorted.vcf.gz",
            toolname=config["behaviors"]["snv-caller"],
        ),
    output:
        "results/export/{projectid}/{sampleid}_{lsid}_{sqid}.snv.vcf.gz",
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
        'bcftools annotate -h <(echo "##wgs-pipelineVersion={params.pipeline_version}") -O u {input} | '
        'bcftools view -i \'(FILTER = "PASS" | FILTER = ".")\' -O u | '
        "bcftools norm -m -both -O u | "
        "bcftools view -i 'FORMAT/DP >= 10 & FORMAT/GQ >= 20 & "
        '((FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 0.2 & FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.8 & GT != "1/1") | '
        ' (FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.05 & GT = "1/1"))\' -O v | '
        'bcftools reheader -s <(echo -e "{wildcards.sqid}\\t{params.exportid}") | bgzip -c > {output}'


rule checksum:
    """
    Create checksum files to go along with transported files
    """
    input:
        "{prefix}",
    output:
        "{prefix}.md5",
    threads: 1
    resources:
        mem_mb="1000",
        qname="small",
    shell:
        "md5sum {input} > {output}"


localrules:
    create_export_manifest,


rule create_export_manifest:
    """
    For export tracking and checkpointing purposes,
    report the expected fileset for a data release
    """
    input:
        bam=lambda wildcards: construct_export_files(wildcards, manifest, "bam"),
        bai=lambda wildcards: construct_export_files(wildcards, manifest, "bai"),
        vcf=lambda wildcards: construct_export_files(wildcards, manifest, "snv.vcf.gz"),
        tbi=lambda wildcards: construct_export_files(wildcards, manifest, "snv.vcf.gz.tbi"),
        bam_md5=lambda wildcards: construct_export_files(wildcards, manifest, "bam.md5"),
        bai_md5=lambda wildcards: construct_export_files(wildcards, manifest, "bai.md5"),
        vcf_md5=lambda wildcards: construct_export_files(wildcards, manifest, "snv.vcf.gz.md5"),
        tbi_md5=lambda wildcards: construct_export_files(wildcards, manifest, "snv.vcf.gz.tbi.md5"),
    output:
        "results/export/{projectid}/manifest.tsv",
    shell:
        "echo {input.bam} {input.bai} {input.vcf} {input.tbi} | sed 's/ /\\n/g' > {output}"
