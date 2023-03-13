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
    container:
        "docker://rocker/tidyverse:latest"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    script:
        "../scripts/construct_linker_from_labbook.R"


rule create_cram_export:
    """
    Take bqsr bamfile and turn it into something to release
    """
    input:
        bam="results/bqsr/{projectid}/{sqid}.bam",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        "results/export/{projectid}/{sampleid}_{lsid}_{sqid}.cram",
    params:
        pipeline_version=pipeline_version,
        reference_url=config["references"][reference_build]["fasta"],
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_cram_export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/samtools.yaml"
    container:
        "{}/bwa.sif".format(apptainer_images)
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "samtools reheader -c 'sed \"s/SM:{wildcards.sqid}/SM:{params.exportid}/ ; "
        "s/LB:{wildcards.sqid}/LB:{params.exportid}/ ; "
        "s/PU:{wildcards.sqid}/PU:{params.exportid}/ ; "
        "\\$a@CO\\twgs-pipelineVersion={params.pipeline_version}\\n@CO\\treferenceUrl={params.reference_url}\"' {input.bam} | "
        "samtools view -C -T {input.fasta} -o {output}"


rule create_crai_export:
    """
    Index export cram file

    This isn't using rule inheritance from the rule in the picard code
    due to idiosyncratic requirements of loading order of rules with inheritance
    """
    input:
        cram="results/export/{prefix}.cram",
    output:
        crai="results/export/{prefix}.crai",
    benchmark:
        "results/performance_benchmarks/create_crai_export/{prefix}.tsv"
    conda:
        "../envs/samtools.yaml"
    container:
        "{}/bwa.sif".format(apptainer_images)
    threads: 4
    resources:
        mem_mb="8000",
        qname="small",
    shell:
        "samtools index -@ {threads} -o {output.crai} {input.cram}"


rule create_snv_vcf_export:
    """
    Take snv vcf output and turn it into something to release

    Note the following things that have nothing to do with actual vcf spec compliance:

    - Moon, the first intended downstream user of this file, has some truly absurd logic
    for determining reference genome build that involves sniffing the vcf header for the
    case-sensitive strings 'GRCh38' 'hg38' 'GRCh37' 'hg19'. as such, the user configuration
    genome build tag is modified to try to meet that format restriction

    - Moon commits the cardinal sin of reimplementing a vcf parser. It does not seem
    to bother to implement a handler for unphased heterozygotes encoded as '1/0', even
    though they are completely valid and equivalent to hets encoded as '0/1'. I don't
    know at this time whether this is actually causing any particular issue with Moon,
    but I'm going to resentfully change these genotypes manually in anticipation.
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sqid}}.sorted.vcf.gz",
            toolname=config["behaviors"]["snv-caller"],
        ),
    output:
        temp("results/export/{projectid}/{sampleid}_{lsid}_{sqid}.snv-allregions.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_snv_vcf_export/export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    container:
        "{}/bcftools.sif".format(apptainer_images)
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        'bcftools annotate -h <(echo -e "##wgs-pipelineVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O u {input} | '
        'bcftools view -i \'(FILTER = "PASS" | FILTER = ".")\' -O u | '
        "bcftools norm -m -both -O u | "
        "bcftools view -i 'FORMAT/DP >= 10 & FORMAT/GQ >= 20 & "
        '((FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 0.2 & FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.8 & GT != "1/1") | '
        ' (FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.05 & GT = "1/1"))\' -O v | '
        'bcftools reheader -s <(echo -e "{wildcards.sqid}\\t{params.exportid}") | '
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


use rule create_snv_vcf_export as create_snv_vcf_nonexport with:
    output:
        temp("results/nonexport/{projectid}/{sqid}.snv-allregions.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sqid}",
    benchmark:
        "results/performance_benchmarks/create_snv_vcf_export/nonexport/{projectid}/{sqid}.tsv"


rule remove_snv_region_exclusions:
    """
    Once SNV output data have had hard filters applied, further remove configurable exclusion regions.
    These are intended to be pulled from https://github.com/Boyle-Lab/Blacklist
    """
    input:
        vcf="{prefix}.snv-allregions.vcf.gz",
        bed="reference_data/references/{}/ref.exclusion.regions.bed".format(reference_build),
    output:
        vcf="{prefix}.snv.vcf.gz",
    benchmark:
        "results/remove_snv_region_exclusions/{prefix}.tsv"
    conda:
        "../envs/bedtools.yaml"
    container:
        "{}/bedtools.sif".format(apptainer_images)
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "bedtools intersect -a {input.vcf} -b {input.bed} -wa -v -header | bgzip -c > {output}"


rule create_sv_vcf_export:
    """
    Take sv vcf output and turn it into something to release

    Note the following things that have nothing to do with actual vcf spec compliance:

    - Moon, the first intended downstream user of this file, has some truly absurd logic
    for determining reference genome build that involves sniffing the vcf header for the
    case-sensitive strings 'GRCh38' 'hg38' 'GRCh37' 'hg19'. as such, the user configuration
    genome build tag is modified to try to meet that format restriction

    - Moon commits the cardinal sin of reimplementing a vcf parser. It does not seem
    to bother to implement a handler for unphased heterozygotes encoded as '1/0', even
    though they are completely valid and equivalent to hets encoded as '0/1'. I don't
    know at this time whether this is actually causing any particular issue with Moon,
    but I'm going to resentfully change these genotypes manually in anticipation.
    """
    input:
        "results/final/{projectid}/{sqid}.sv.vcf.gz",
    output:
        "results/export/{projectid}/{sampleid}_{lsid}_{sqid}.sv.vcf.gz",
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_sv_vcf_export/export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    container:
        "{}/bcftools.sif".format(apptainer_images)
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        'bcftools annotate -h <(echo -e "##wgs-pipelineVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O v {input} | '
        'bcftools reheader -s <(echo -e "{wildcards.sqid}\\t{params.exportid}") | '
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


use rule create_sv_vcf_export as create_sv_vcf_nonexport with:
    output:
        temp("results/nonexport/{projectid}/{sqid}.sv.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sqid}",
    benchmark:
        "results/performance_benchmarks/create_sv_vcf_export/nonexport/{projectid}/{sqid}.tsv"


rule checksum:
    """
    Create checksum files to go along with transported files
    """
    input:
        "results/{export_status}/{projectid}/{prefix}",
    output:
        "results/{export_status}/{projectid}/{prefix}.md5",
    threads: 1
    resources:
        mem_mb="1000",
        qname="small",
    shell:
        "md5sum {input} | sed 's|results/export/{wildcards.projectid}/||' > {output}"


localrules:
    create_export_manifest,
    create_nonexported_manifest,


rule create_export_manifest:
    """
    For export tracking and checkpointing purposes,
    report the expected fileset for a data release
    """
    input:
        cram=lambda wildcards: ed.construct_export_files(wildcards, manifest, checkpoints, "cram"),
        crai=lambda wildcards: ed.construct_export_files(wildcards, manifest, checkpoints, "crai"),
        vcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz"
        ),
        tbi=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz.tbi"
        ),
        sv_vcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz"
        ),
        sv_tbi=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.tbi"
        ),
        cram_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "cram.md5"
        ),
        crai_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "crai.md5"
        ),
        vcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz.md5"
        ),
        tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz.tbi.md5"
        ),
        sv_vcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.md5"
        ),
        sv_tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.tbi.md5"
        ),
    output:
        "results/export/{projectid}/manifest.tsv",
    shell:
        "echo {input.cram} {input.crai} {input.vcf} {input.tbi} {input.sv_vcf} {input.sv_tbi} | sed 's/ /\\n/g' > {output}"


rule create_nonexported_manifest:
    """
    At least one library in each run is expected
    to not be exported. I'm still interested in those
    libraries' variant calling performance. Don't bother
    constructing crams, as they're gigantic and annoying.
    """
    input:
        vcf=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz"
        ),
        tbi=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz.tbi"
        ),
        vcf_md5=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz.md5"
        ),
        tbi_md5=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz.tbi.md5"
        ),
        sv_vcf=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz"
        ),
        sv_tbi=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.tbi"
        ),
        sv_vcf_md5=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.md5"
        ),
        sv_tbi_md5=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.tbi.md5"
        ),
    output:
        "results/nonexport/{projectid}/manifest.tsv",
    shell:
        "echo {input.vcf} {input.tbi} {input.sv_vcf} {input.sv_tbi} | sed 's/ /\\n/g' > {output}"


rule export_data:
    """
    Lossily move results/export data contents to deployment directory
    somewhere else
    """
    input:
        "workflow/scripts/export_data.bash",
    output:
        "results/export/md5_checks.txt",
    params:
        export_directory=config["behaviors"]["export-directory"],
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "{input} {params.export_directory} {output}"
