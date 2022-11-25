rule duphold_run:
    """
    Use duphold to annotate SVs in a vcf with various quality guesses
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        snv_vcf="results/octopus/{projectid}/{sampleid}.sorted.vcf.gz",
        snv_tbi="results/octopus/{projectid}/{sampleid}.sorted.vcf.gz.tbi",
        sv_vcf="results/{toolname}/{projectid}/{sampleid}.{toolname}.vcf.gz",
        sv_tbi="results/{toolname}/{projectid}/{sampleid}.{toolname}.vcf.gz.tbi",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        bcf=temp("results/{toolname}/{projectid}/{sampleid}.duphold-annotated.bcf"),
    conda:
        "../envs/duphold.yaml"
    threads: 4
    resources:
        h_vmem="8000",
        qname="small",
    shell:
        "duphold -s {input.snv_vcf} -t {threads} -v {input.sv_vcf} -b {input.bam} -f {input.fasta} -o {output.bcf}"


rule duphold_apply:
    """
    After running duphold, actually use the filters to get a (hopefully) cleaner dataset
    """
    input:
        bcf="results/{toolname}/{projectid}/{sampleid}.duphold-annotated.bcf",
    output:
        vcf="results/{toolname}/{projectid}/{sampleid}.duphold-filtered.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        h_vmem="4000",
        qname="small",
    shell:
        'bcftools view -i \'(SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] > 1.3)\' '
        "--threads {threads} -O z -o {output.vcf} {input.bcf}"
