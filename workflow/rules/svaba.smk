rule svaba_run:
    """
    Run svaba in germline mode
    """
    input:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
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
        "../envs/svaba.yaml"
    threads: 8
    resources:
        mem_mb="32000",
        qname="large",
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
    """
    input:
        vcf="results/svaba/{projectid}/{sampleid}.svaba.unfiltered.sv.vcf",
    output:
        vcf="results/svaba/{projectid}/{sampleid}.svaba.vcf.gz",
        linker=temp("results/svaba/{projectid}/{sampleid}.svaba.reheader_linker.tsv"),
    params:
        tmpdir="temp",
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        blacklist_definition='##FILTER=<ID=BLACKLIST,Description=\\"Variant is in calling exclusion region\\">',
    benchmark:
        "results/performance_benchmarks/svaba_select_output_variants/{projectid}/{sampleid}.svaba.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb="16000",
        qname="small",
    shell:
        "mkdir -p {params.tmpdir} && "
        "echo -e '{params.bam}\\t{wildcards.sampleid}' > {output.linker} && "
        'awk -v blacklist="{params.blacklist_definition}" \'/^#CHROM/ {{OFS="\\t" ; print blacklist"\\n"$0}} ; ! /^#CHROM/\' {input.vcf} | '
        'awk -F"\\t" \'/^#/ ; ! /^#/ {{OFS = "\\t" ; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13}}\' | '
        "sed 's/<ID=PL,Number=.,/<ID=PL,Number=G,/' | "
        "bcftools reheader -s {output.linker} | "
        "bcftools sort -O z --temp-dir {params.tmpdir} -o {output.vcf}"
