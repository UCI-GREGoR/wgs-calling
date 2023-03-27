rule duphold_run:
    """
    Use duphold to annotate SVs in a vcf with various quality guesses
    """
    input:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
        snv_vcf=expand(
            "results/{caller}/{{projectid}}/{{sampleid}}.sorted.vcf.gz",
            caller=config["behaviors"]["snv-caller"],
        ),
        snv_tbi=expand(
            "results/{caller}/{{projectid}}/{{sampleid}}.sorted.vcf.gz.tbi",
            caller=config["behaviors"]["snv-caller"],
        ),
        sv_vcf="results/{toolname}/{projectid}/{sampleid}.{toolname}.vcf.gz",
        sv_tbi="results/{toolname}/{projectid}/{sampleid}.{toolname}.vcf.gz.tbi",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
    output:
        bcf=temp("results/{toolname}/{projectid}/{sampleid}.{toolname}.duphold-annotated.bcf"),
    benchmark:
        "results/performance_benchmarks/duphold_run/{toolname}/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/duphold.yaml"
    threads: config_resources["duphold"]["threads"]
    resources:
        mem_mb=config_resources["duphold"]["memory"],
        qname=rc.select_queue(config_resources["duphold"]["queue"]),
    shell:
        "duphold -s {input.snv_vcf} -t {threads} -v {input.sv_vcf} -b {input.bam} -f {input.fasta} -o {output.bcf}"


rule duphold_apply:
    """
    After running duphold, actually use the filters to get a (hopefully) cleaner dataset
    """
    input:
        bcf="results/{toolname}/{projectid}/{sampleid}.{toolname}.duphold-annotated.bcf",
    output:
        vcf="results/{toolname}/{projectid}/{sampleid}.{toolname}.duphold-filtered.vcf.gz",
    benchmark:
        "results/performance_benchmarks/duphold_apply/{toolname}/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=rc.select_queue(config_resources["bcftools"]["queue"]),
    shell:
        'bcftools view -i \'(FILTER = "PASS" | FILTER = ".") & '
        '((FMT/DHFFC[0] = ".") | '
        ' (SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | '
        ' (SVTYPE != "DEL" & SVTYPE != "INS" & FMT/DHBFC[0] > 1.3) | '
        ' (SVTYPE = "INS")) \' '
        "--threads {threads} -O z -o {output.vcf} {input.bcf}"
