rule truvari_merge_within_caller:
    """
    Use truvari to collapse approximately redundant variants within a single caller's output
    """
    input:
        vcf="results/{toolname}/{projectid}/{sampleid}.{toolname}.duphold-filtered.vcf.gz",
        tbi="results/{toolname}/{projectid}/{sampleid}.{toolname}.duphold-filtered.vcf.gz.tbi",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
    output:
        vcf=temp("results/final/{projectid}/{sampleid}.{toolname}.within-merge.vcf.gz"),
        collapsed=temp(
            "results/final/{projectid}/{sampleid}.{toolname}.within-merge.collapsed.vcf.gz"
        ),
    benchmark:
        "results/performance_benchmarks/truvari_merge_within_caller/{projectid}/{sampleid}.{toolname}.tsv"
    conda:
        "../envs/truvari.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "truvari collapse -i {input.vcf} -o {output.vcf} -c {output.collapsed} -f {input.fasta} "
        "-p 0.5 -O 0.25 -P 0.5"


rule bcftools_concat_sv_callers:
    """
    Use bcftools to concatenate variants from multiple SV callers into a single file,
    with the intention of providing this file to a second round of truvari.
    """
    input:
        vcf=expand(
            "results/final/{{projectid}}/{{sampleid}}.{toolname}.within-merge.vcf.gz",
            toolname=config["behaviors"]["sv-callers"],
        ),
        tbi=expand(
            "results/final/{{projectid}}/{{sampleid}}.{toolname}.within-merge.vcf.gz.tbi",
            toolname=config["behaviors"]["sv-callers"],
        ),
    output:
        vcf=temp("results/final/{projectid}/{sampleid}.tool-sorted.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/bcftools_concat_sv_callers/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "bcftools concat -a -d exact -O u | bcftools sort -O z -o {output.vcf}"


use rule truvari_merge_within_caller as truvari_merge_between_callers with:
    input:
        vcf="results/final/{projectid}/{sampleid}.tool-sorted.vcf.gz",
        tbi="results/final/{projectid}/{sampleid}.tool-sorted.vcf.gz.tbi",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
    output:
        vcf=temp("results/final/{projectid}/{sampleid}.truvari-raw.vcf.gz"),
        collapsed=temp("results/final/{projectid}/{sampleid}.truvari-raw.collapsed.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/truvari_merge_between_callers/{projectid}/{sampleid}.tsv"


rule truvari_ensemble_sv_vcf:
    """
    Given a truvari-merged version of SV calls, select variants
    that are flagged as present in at least some user-configurable
    number of SV caller inputs
    """
    input:
        "results/final/{projectid}/{sampleid}.sv.truvari-raw.vcf.gz",
    output:
        "results/final/{projectid}/{sampleid}.sv.vcf.gz",
    params:
        bcftools_filter_count="INFO/FOUNDBY = {}".format(
            config["behaviors"]["sv-ensemble"]["min-count"]
        ),
        bcftools_filter_sources=" & INFO/svdb_origin ~ '"
        + "' & INFO/svdb_origin ~ '".join(config["behaviors"]["sv-ensemble"]["required-callers"])
        + "'"
        if "required-callers" in config["behaviors"]["sv-ensemble"]
        else "",
    benchmark:
        "results/performance_benchmarks/truvari_ensemble_sv_vcf/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        'bcftools filter -i "{params.bcftools_filter_count} {params.bcftools_filter_sources}" -O z -o {output} {input}'


rule truvari_summarize_sv_variant_sources:
    """
    For a raw merge output from truvari, use bcftools to extract reasonable summary information about
    the variant from its info content, with the goal of getting this information into an Rmd report
    """
    input:
        "results/final/{projectid}/{sampleid}.sv.truvari-raw.vcf.gz",
    output:
        "results/reports/sv_data/{projectid}/{sampleid}.sv.truvari-raw.tsv",
    benchmark:
        "results/performance_benchmarks/truvari_summarize_sv_variant_sources/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/svdb_origin\\n' {input} > {output}"
