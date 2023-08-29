checkpoint generate_linker:
    """
    From a sample logbook, generate a simple linker
    between various sample ID types
    """
    output:
        linker="results/export/linker.tsv",
    params:
        logbook=config["sample-logbook"] if "sample-logbook" in config else None,
        sex_linker=config["sample-linking"]["sex"]
        if "sample-linking" in config and "sex" in config["sample-linking"]
        else None,
        external_id_linker=config["sample-linking"]["external-ids"]
        if "sample-linking" in config and "external-ids" in config["sample-linking"]
        else None,
    benchmark:
        "results/performance_benchmarks/generate_linker/linker.tsv"
    conda:
        "../envs/r.yaml" if not use_containers else None
    container:
        "{}/r.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["r"]["threads"]
    resources:
        mem_mb=config_resources["r"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["r"]["queue"], config_resources["queues"]
        ),
    script:
        "../scripts/construct_linker_from_inputs.R"


rule create_cram_export:
    """
    Take (conceptually) post-bqsr bamfile and turn it into something to release
    """
    input:
        bam="results/aligned_bams/{projectid}/{sqid}.bam",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        temp("results/export/{projectid}/{sampleid}_{lsid}_{sqid}.cram"),
    params:
        pipeline_version=pipeline_version,
        reference_url=config["references"][reference_build]["fasta"],
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_cram_export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/bwa.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools"]["threads"]
    resources:
        mem_mb=config_resources["samtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "samtools reheader -c 'sed \"s/SM:{wildcards.sqid}/SM:{params.exportid}/ ; "
        "s/LB:{wildcards.sqid}/LB:{params.exportid}/ ; "
        "s/PU:{wildcards.sqid}/PU:{params.exportid}/ ; "
        "\\$a@CO\\twgs-pipelineVersion={params.pipeline_version}\\n@CO\\treferenceUrl={params.reference_url}\"' {input.bam} | "
        "samtools view -C -T {input.fasta} -o {output}"


rule create_cram_export_simplified_id:
    """
    Take (conceptually) post-bqsr bamfile and turn it into something to release.
    The prior version of this rule was designed for use with old-style
    labbook annotations including analyte ID, etc. This streamlined version
    removes some of the distinctiveness of IDs in favor of streamlined output.
    """
    input:
        bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        temp("results/export/{projectid}/{sampleid}.cram"),
    params:
        pipeline_version=pipeline_version,
        reference_url=config["references"][reference_build]["fasta"],
    benchmark:
        "results/performance_benchmarks/create_cram_export_simplified_id/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/bwa.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools"]["threads"]
    resources:
        mem_mb=config_resources["samtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "samtools reheader -c 'sed -r \""
        "s/SM:[^ \\t]+([ \\t]+)/SM:{wildcards.sampleid}\\1/ ; "
        "s/LB:[^ \\t]+([ \\t]+)/LB:{wildcards.sampleid}\\1/ ; "
        "s/PU:[^ \\t]+([ \\t]+)/PU:{wildcards.sampleid}\\1/ ; "
        "s/SM:[^ \\t]+$/SM:{wildcards.sampleid}/ ; "
        "s/LB:[^ \\t]+$/LB:{wildcards.sampleid}/ ; "
        "s/PU:[^ \\t]+$/PU:{wildcards.sampleid}/ ; "
        "\\$a@CO\\twgs-pipelineVersion={params.pipeline_version}\\n@CO\\treferenceUrl={params.reference_url}\"' {input.bam} | "
        "samtools view -C -T {input.fasta} -o {output}"


rule create_crai_export:
    """
    Index export cram file.

    This isn't using rule inheritance from the rule in the picard code
    due to idiosyncratic requirements of loading order of rules with inheritance.
    """
    input:
        cram="results/export/{prefix}.cram",
    output:
        crai=temp("results/export/{prefix}.crai"),
    benchmark:
        "results/performance_benchmarks/create_crai_export/{prefix}.tsv"
    conda:
        "../envs/samtools.yaml" if not use_containers else None
    container:
        "{}/bwa.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["samtools"]["threads"]
    resources:
        mem_mb=config_resources["samtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "samtools index -@ {threads} -o {output.crai} {input.cram}"


rule create_snv_gvcf_export:
    """
    Take snv g.vcf output and turn it into something to release.

    *Some* of the modifications applied to snv vcfs are applied here, but
    filtering is deferred to downstream calling.
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sqid}}.sorted.g.vcf.gz",
            toolname=config["behaviors"]["snv-caller"],
        ),
    output:
        temp("results/export/{projectid}/{sampleid}_{lsid}_{sqid}.snv.g.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_snv_gvcf_export/export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        'bcftools annotate -h <(echo -e "##wgs-pipelineVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O v {input} | '
        'bcftools reheader -s <(echo -e "{wildcards.sqid}\\t{params.exportid}") | '
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


use rule create_snv_gvcf_export as create_snv_gvcf_nonexport with:
    output:
        "results/nonexport/{projectid}/{sqid}.snv.g.vcf.gz",
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sqid}",
    benchmark:
        "results/performance_benchmarks/create_snv_vcf_export/nonexport/{projectid}/{sqid}.tsv"


rule create_snv_gvcf_export_simplified_id:
    """
    Take snv g.vcf output and turn it into something to release.

    *Some* of the modifications applied to snv vcfs are applied here, but
    filtering is deferred to downstream calling.
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sampleid}}.sorted.g.vcf.gz",
            toolname=config["behaviors"]["snv-caller"],
        ),
    output:
        temp("results/export/{projectid}/{sampleid}.snv.g.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}",
    benchmark:
        "results/performance_benchmarks/create_snv_gvcf_export_simplified_id/export/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        'bcftools annotate -h <(echo -e "##wgs-pipelineVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O v {input} | '
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


use rule create_snv_gvcf_export_simplified_id as create_snv_gvcf_nonexport_simplified_id with:
    output:
        "results/nonexport/{projectid}/{sampleid}.snv.g.vcf.gz",
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}",
    benchmark:
        "results/performance_benchmarks/create_snv_vcf_export_simplified_id/nonexport/{projectid}/{sampleid}.tsv"


rule create_snv_vcf_export:
    """
    Take snv vcf output and turn it into something to release.

    Note the following things that have nothing to do with actual vcf spec compliance:

    - for some tool compatibility, modify genome build code to standardized format

    - some non-compliant vcf parsers cannot handle the valid unphased heterozygote '1/0',
    and as such the genotypes are "simplified" to '0/1'
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
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
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


rule create_snv_vcf_export_simplified_id:
    """
    Take snv vcf output and turn it into something to release, but with a slightly
    simpler ID format.

    Note the following things that have nothing to do with actual vcf spec compliance:

    - for some tool compatibility, modify genome build code to standardized format

    - some non-compliant vcf parsers cannot handle the valid unphased heterozygote '1/0',
    and as such the genotypes are "simplified" to '0/1'
    """
    input:
        expand(
            "results/{toolname}/{{projectid}}/{{sampleid}}.sorted.vcf.gz",
            toolname=config["behaviors"]["snv-caller"],
        ),
    output:
        temp("results/export/{projectid}/{sampleid}.snv-allregions.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}",
    benchmark:
        "results/performance_benchmarks/create_snv_vcf_export_simplified_id/export/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        'bcftools annotate -h <(echo -e "##wgs-pipelineVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O u {input} | '
        'bcftools view -i \'(FILTER = "PASS" | FILTER = ".")\' -O u | '
        "bcftools norm -m -both -O u | "
        "bcftools view -i 'FORMAT/DP >= 10 & FORMAT/GQ >= 20 & "
        '((FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 0.2 & FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.8 & GT != "1/1") | '
        ' (FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.05 & GT = "1/1"))\' -O v | '
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


use rule create_snv_vcf_export_simplified_id as create_snv_vcf_nonexport_simplified_id with:
    output:
        temp("results/nonexport/{projectid}/{sampleid}.snv-allregions.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}",
    benchmark:
        "results/performance_benchmarks/create_snv_vcf_export_simplified_id/nonexport/{projectid}/{sampleid}.tsv"


rule remove_snv_region_exclusions_export:
    """
    Once SNV output data have had hard filters applied, further remove configurable exclusion regions.
    These are intended to be pulled from https://github.com/Boyle-Lab/Blacklist
    """
    input:
        vcf="results/export/{projectid}/{sqid}.snv-allregions.vcf.gz",
        bed="reference_data/references/{}/ref.exclusion.regions.bed".format(reference_build),
    output:
        vcf=temp("results/export/{projectid}/{sqid}.snv.vcf.gz"),
    benchmark:
        "results/remove_snv_region_exclusions_export/{projectid}/{sqid}.tsv"
    conda:
        "../envs/bedtools.yaml" if not use_containers else None
    container:
        "{}/bedtools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bedtools"]["threads"]
    resources:
        mem_mb=config_resources["bedtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bedtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bedtools intersect -a {input.vcf} -b {input.bed} -wa -v -header | bgzip -c > {output}"


use rule remove_snv_region_exclusions_export as remove_snv_region_exclusions_nonexport with:
    input:
        vcf="results/nonexport/{projectid}/{sqid}.snv-allregions.vcf.gz",
        bed="reference_data/references/{}/ref.exclusion.regions.bed".format(reference_build),
    output:
        vcf="results/nonexport/{projectid}/{sqid}.snv.vcf.gz",
    benchmark:
        "results/remove_snv_region_exclusions_nonexport/{projectid}/{sqid}.tsv"


rule create_sv_vcf_export:
    """
    Take sv vcf output and turn it into something to release

    Note the following things that have nothing to do with actual vcf spec compliance:

    - for some tool compatibility, modify genome build code to standardized format

    - some non-compliant vcf parsers cannot handle the valid unphased heterozygote '1/0',
    and as such the genotypes are "simplified" to '0/1'
    """
    input:
        "results/final/{projectid}/{sqid}.sv.vcf.gz",
    output:
        temp("results/export/{projectid}/{sampleid}_{lsid}_{sqid}.sv.with-bnd.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}_{lsid}_{sqid}",
    benchmark:
        "results/performance_benchmarks/create_sv_vcf_export/export/{projectid}/{sampleid}_{lsid}_{sqid}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        'bcftools annotate -h <(echo -e "##wgs-pipelineVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O v {input} | '
        'bcftools reheader -s <(echo -e "{wildcards.sqid}\\t{params.exportid}") | '
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


use rule create_sv_vcf_export as create_sv_vcf_nonexport with:
    output:
        "results/nonexport/{projectid}/{sqid}.sv.with-bnd.vcf.gz",
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sqid}",
    benchmark:
        "results/performance_benchmarks/create_sv_vcf_export/nonexport/{projectid}/{sqid}.tsv"


rule create_sv_vcf_export_simplified_id:
    """
    Take sv vcf output and turn it into something to release, but with a slightly
    simpler ID format

    Note the following things that have nothing to do with actual vcf spec compliance:

    - for some tool compatibility, modify genome build code to standardized format

    - some non-compliant vcf parsers cannot handle the valid unphased heterozygote '1/0',
    and as such the genotypes are "simplified" to '0/1'
    """
    input:
        "results/final/{projectid}/{sampleid}.sv.vcf.gz",
    output:
        temp("results/export/{projectid}/{sampleid}.sv.with-bnd.vcf.gz"),
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}",
    benchmark:
        "results/performance_benchmarks/create_sv_vcf_export_simplified_id/export/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        'bcftools annotate -h <(echo -e "##wgs-pipelineVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O v {input} | '
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


use rule create_sv_vcf_export_simplified_id as create_sv_vcf_nonexport_simplified_id with:
    output:
        "results/nonexport/{projectid}/{sampleid}.sv.with-bnd.vcf.gz",
    params:
        pipeline_version=pipeline_version,
        reference_build=lambda wildcards: sm.format_reference_build(reference_build),
        exportid="{sampleid}",
    benchmark:
        "results/performance_benchmarks/create_sv_vcf_export/nonexport/{projectid}/{sampleid}.tsv"


rule remove_breakends:
    """
    Conditionally remove breakends from ensemble called SVs based on user configuration
    """
    input:
        "{prefix}.sv.with-bnd.vcf.gz",
    output:
        temp("{prefix}.sv.vcf.gz"),
    params:
        remove_breakends=config["behaviors"]["sv-remove-breakends"],
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        'if [[ "{params.remove_breakends}" == "True" ]] ; then '
        "bcftools filter -i 'SVTYPE != \"BND\"' -O z -o {output} {input} ; else "
        "cp {input} {output} ; "
        "fi"


rule checksum_export:
    """
    Create checksum files to go along with transported files
    """
    input:
        "results/export/{projectid}/{prefix}",
    output:
        temp("results/export/{projectid}/{prefix}.md5"),
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "md5sum {input} | sed -r 's|  .*/([^/ ]+)$|  \\1|' > {output}"


use rule checksum_export as checksum_nonexport with:
    input:
        "results/nonexport/{projectid}/{prefix}",
    output:
        "results/nonexport/{projectid}/{prefix}.md5",


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
        gvcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz"
        ),
        gvcf_tbi=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi"
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
        gvcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.md5"
        ),
        gvcf_tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi.md5"
        ),
        sv_vcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.md5"
        ),
        sv_tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.tbi.md5"
        ),
    output:
        temp("results/export/{projectid}/manifest.tsv"),
    shell:
        "echo {input.cram} {input.crai} {input.vcf} {input.tbi} "
        "{input.gvcf} {input.gvcf_tbi} {input.sv_vcf} {input.sv_tbi} | "
        "sed 's/ /\\n/g' > {output}"


rule create_nonexported_manifest:
    """
    Periodically, one or more libraries in a run are
    not expected to be exported. Regardless, we might still
    be interested in those libraries' variant calling
    performance. Don't bother constructing crams, as they're
    gigantic and annoying.
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
        gvcf=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz"
        ),
        gvcf_tbi=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi"
        ),
        gvcf_md5=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.md5"
        ),
        gvcf_tbi_md5=lambda wildcards: ed.construct_nonexport_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi.md5"
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
        temp("results/nonexport/{projectid}/manifest.tsv"),
    shell:
        "echo {input.vcf} {input.tbi} {input.gvcf} {input.gvcf_tbi} "
        "{input.sv_vcf} {input.sv_tbi} | sed 's/ /\\n/g' > {output}"


rule zip_vcfs:
    """
    At downstream user request, zip all vcfs without indices for single download.
    """
    input:
        vcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz"
        ),
        sv_vcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz"
        ),
    output:
        "results/export/{projectid}/{projectid}_vcfs.zip",
    shell:
        "zip -j {output} {input.vcf} {input.sv_vcf}"


rule export_data_local:
    """
    Move results/export data contents to deployment directory
    somewhere else
    """
    input:
        bash="workflow/scripts/export_data.bash",
        cram=lambda wildcards: ed.construct_export_files(wildcards, manifest, checkpoints, "cram"),
        crai=lambda wildcards: ed.construct_export_files(wildcards, manifest, checkpoints, "crai"),
        vcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz"
        ),
        tbi=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.vcf.gz.tbi"
        ),
        gvcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz"
        ),
        gvcf_tbi=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi"
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
        gvcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.md5"
        ),
        gvcf_tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi.md5"
        ),
        sv_vcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.md5"
        ),
        sv_tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.tbi.md5"
        ),
        vcf_zip=expand(
            "results/export/{projectid}/{projectid}_vcfs.{suffix}", suffix=["zip", "zip.md5"]
        ),
        manifest="results/export/{projectid}/manifest.tsv",
    output:
        "results/export/{projectid}/md5_checks.txt",
    params:
        export_directory=config["behaviors"]["export-directory"]
        if "export-directory" in config["behaviors"]
        else "",
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "{input.bash} {params.export_directory} {output}"


rule export_data_remote:
    """
    Sync results/export data contents to remote deployment s3 bucket
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
        gvcf=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz"
        ),
        gvcf_tbi=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi"
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
        gvcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.md5"
        ),
        gvcf_tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "snv.g.vcf.gz.tbi.md5"
        ),
        sv_vcf_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.md5"
        ),
        sv_tbi_md5=lambda wildcards: ed.construct_export_files(
            wildcards, manifest, checkpoints, "sv.vcf.gz.tbi.md5"
        ),
        manifest="results/export/{projectid}/manifest.tsv",
    output:
        "results/export/{projectid}/s3_transfer_complete.txt",
    params:
        export_dir="results/export/{projectid}",
        bucketname=config["behaviors"]["export-s3"]["bucket-name"]
        if "export-s3" in config["behaviors"]
        else None,
        profile="--profile {}".format(config["behaviors"]["export-s3"]["profile-name"])
        if "export-s3" in config["behaviors"]
        and "profile-name" in config["behaviors"]["export-s3"]
        else "",
    threads: config_resources["awscli"]["threads"]
    resources:
        mem_mb=config_resources["awscli"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["awscli"]["queue"], config_resources["queues"]
        ),
    shell:
        'aws s3 sync {params.profile} --exclude="*" --include="*.cram*" --include="*.crai*" {params.export_dir} {params.bucketname}/wgs-short-read/{wildcards.projectid}/crams && '
        'aws s3 sync {params.profile} --exclude="*" --include="*.snv.vcf*" {params.export_dir} {params.bucketname}/wgs-short-read/{wildcards.projectid}/snv_vcfs && '
        'aws s3 sync {params.profile} --exclude="*" --include="*.snv.g.vcf*" {params.export_dir} {params.bucketname}/wgs-short-read/{wildcards.projectid}/snv_gvcfs && '
        'aws s3 sync {params.profile} --exclude="*" --include="*.sv.vcf*" {params.export_dir} {params.bucketname}/wgs-short-read/{wildcards.projectid}/sv_vcfs && '
        "touch {output}"


rule export_fastqs_remote:
    """
    Sync fastqs to remote deployment s3 bucket

    I'm not currently certain that this is exactly how I want this process to work,
    so consider this rule somewhat WIP
    """
    input:
        fastqs_r1=lambda wildcards: [
            "results/fastqs/{}/{}_{}_R1_001.fastq.gz".format(wildcards.projectid, sampleid, lane)
            for sampleid, lane in zip(
                manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"],
                manifest.loc[manifest["projectid"] == wildcards.projectid, "lane"],
            )
        ],
        fastqs_r2=lambda wildcards: [
            "results/fastqs/{}/{}_{}_R2_001.fastq.gz".format(wildcards.projectid, sampleid, lane)
            for sampleid, lane in zip(
                manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"],
                manifest.loc[manifest["projectid"] == wildcards.projectid, "lane"],
            )
        ],
    output:
        "results/fastqs/{projectid}/s3_transfer_complete.txt",
    params:
        export_dir="results/fastqs/{projectid}",
        bucketname=config["behaviors"]["export-s3"]["bucket-name"]
        if "export-s3" in config["behaviors"]
        else None,
        profile="--profile {}".format(config["behaviors"]["export-s3"]["profile-name"])
        if "export-s3" in config["behaviors"]
        and "profile-name" in config["behaviors"]["export-s3"]
        else "",
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        'aws s3 sync {params.profile} --exclude="*" --include="*.fastq.gz" {params.export_dir} {params.bucketname}/wgs-short-read/{wildcards.projectid}/fastqs && '
        "touch {output}"
