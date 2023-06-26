rule svaba_run:
    """
    Run svaba in germline mode
    """
    input:
        bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        bai="results/aligned_bams/{projectid}/{sampleid}.bai",
        bwa_fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        bwa_other_files=expand(
            "reference_data/{aligner}/{genome}/ref.fasta.{suffix}",
            aligner=config["behaviors"]["aligner"],
            genome=reference_build,
            suffix=aligner_index_suffixes[config["behaviors"]["aligner"]] + ["fai", "dict"],
        ),
        bed="reference_data/svaba/{}/ref.exclude.bed".format(reference_build),
    output:
        bps="results/svaba/{projectid}/{sampleid}.bps.txt.gz",
        contigs="results/svaba/{projectid}/{sampleid}.contigs.bam",
        discordants="results/svaba/{projectid}/{sampleid}.discordant.txt.gz",
        log="results/svaba/{projectid}/{sampleid}.log",
        alignments="results/svaba/{projectid}/{sampleid}.alignments.txt.gz",
        vcf_indel_unfiltered="results/svaba/{projectid}/{sampleid}.svaba.unfiltered.indel.vcf",
        vcf_sv_unfiltered="results/svaba/{projectid}/{sampleid}.svaba.unfiltered.sv.vcf",
        vcf_indel_filtered="results/svaba/{projectid}/{sampleid}.svaba.indel.vcf",
        vcf_sv_filtered="results/svaba/{projectid}/{sampleid}.svaba.sv.vcf",
    benchmark:
        "results/performance_benchmarks/svaba_run/{projectid}/{sampleid}.svaba.tsv"
    params:
        outprefix="results/svaba/{projectid}/{sampleid}",
    conda:
        "../envs/svaba.yaml" if not use_containers else None
    container:
        "{}/svaba.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["svaba"]["threads"]
    resources:
        mem_mb=config_resources["svaba"]["memory"],
        qname=rc.select_queue(config_resources["svaba"]["queue"], config_resources["queues"]),
    shell:
        "svaba run -p {threads} -G {input.bwa_fasta} -I -L 6 -t {input.bam} -B {input.bed} -a {params.outprefix}"


rule svaba_select_output_variants:
    """
    svaba emits filtered and unfiltered variant sets, split by svs and indels.
    this will probably have to wait until there's actual output files to be
    certain how this should be handled.

    This command has been complicated by malformations in the vcf output format from SvABA.
    The awk and sed and reheader corrections below are as follows:

    - SvABA emits sample IDs as the path to the bamfile that was provided to it. this is obviously
    a problem as it doesn't respect the actual expected sample ID
    - the program flags variants as FILTER=BLACKLIST, but SvABA does not bother to emit a
    FILTER=BLACKLIST definition in the vcf header
    - the output vcf contains mysterious additional columns between FORMAT and the genotypes.
    this has evidently been patched in a non-release version of SvABA. this makes me sad. I've checked
    the commit in question, and evidently it is safe to just ignore the offending columns, as long as
    that happens before bcftools gets to the file
    - the output vcf FORMAT definition for PL is '.' but bcftools complains that it should be 'G';
    that is, one value per genotype

    Furthermore, with the intention of providing this file to svtools for parsing, entries with
    missing breakend pair seem to cause svtools much distress, and as such, those are removed
    before being passed along.
    """
    input:
        vcf="results/svaba/{projectid}/{sampleid}.svaba.unfiltered.sv.vcf",
    output:
        vcf="results/svaba/{projectid}/{sampleid}.svaba.as_bnd.vcf.gz",
        linker=temp("results/svaba/{projectid}/{sampleid}.svaba.reheader_linker.tsv"),
    params:
        tmpdir=tempDir,
        bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        blacklist_definition='##FILTER=<ID=BLACKLIST,Description=\\"Variant is in calling exclusion region\\">',
    benchmark:
        "results/performance_benchmarks/svaba_select_output_variants/{projectid}/{sampleid}.svaba.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=rc.select_queue(config_resources["bcftools"]["queue"], config_resources["queues"]),
    shell:
        "mkdir -p {params.tmpdir} && "
        "echo -e '{params.bam}\\t{wildcards.sampleid}' > {output.linker} && "
        'awk -v blacklist="{params.blacklist_definition}" \'/^#CHROM/ {{OFS="\\t" ; print blacklist"\\n"$0}} ; ! /^#CHROM/\' {input.vcf} | '
        'awk -F"\\t" \'/^#/ ; ! /^#/ {{OFS = "\\t" ; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13}}\' | '
        "sed 's/<ID=PL,Number=.,/<ID=PL,Number=G,/' | "
        "bcftools view -f 'PASS' | "
        "bcftools reheader -s {output.linker} | "
        "bcftools sort -O z --temp-dir {params.tmpdir} -o {output.vcf}"


rule vcf_to_bedpe:
    """
    Use svtools to convert vcf to flattened bedpe format
    """
    input:
        "{prefix}.svaba.as_bnd.vcf.gz",
    output:
        "{prefix}.svaba.as_bnd.bedpe",
    conda:
        "../envs/svtools.yaml" if not use_containers else None
    container:
        "docker://halllab/svtools:v0.5.1" if use_containers else None
    threads: config_resources["svtools"]["threads"]
    resources:
        mem_mb=config_resources["svtools"]["memory"],
        qname=rc.select_queue(config_resources["svtools"]["queue"], config_resources["queues"]),
    shell:
        "gunzip -c {input} > {output}.tmp && "
        "svtools vcftobedpe -i {output}.tmp -o {output} && "
        "rm {output}.tmp"


rule svaba_resolve_breakends:
    """
    Try to represent SvABA BND format variants as other, more specific SV types
    """
    input:
        "{prefix}.svaba.as_bnd.bedpe",
    output:
        "{prefix}.svaba.resolved.bedpe",
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=rc.select_queue(config_resources["default"]["queue"], config_resources["queues"]),
    script:
        "../scripts/reclassify_svs.py"


rule bedpe_to_vcf:
    """
    Use svtools to convert flattened bedpe to vcf format
    """
    input:
        "{prefix}.svaba.resolved.bedpe",
    output:
        "{prefix}.svaba.vcf.gz",
    params:
        tmpdir=tempDir,
    conda:
        "../envs/svtools.yaml" if not use_containers else None
    container:
        "docker://halllab/svtools:v0.5.1" if use_containers else None
    threads: config_resources["svtools"]["threads"]
    resources:
        mem_mb=config_resources["svtools"]["memory"],
        qname=rc.select_queue(config_resources["svtools"]["queue"], config_resources["queues"]),
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        "svtools bedpetovcf -i {input} -o {output}.tmp && "
        "bcftools sort -O z --temp-dir {params.tmpdir} -o {output} {output}.tmp && "
        "rm {output}.tmp"
