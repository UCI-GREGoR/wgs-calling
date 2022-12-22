rule svaba_run:
    """
    Run svaba in germline mode
    """
    input:
        bam="results/markdups/{projectid}/{sampleid}.mrkdup.sort.bam",
        bai="results/bwa-mem2/{projectid}/{sampleid}.mrkdup.sort.bam.bai",
        bwa_fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        bwa_other_files=expand(
            "reference_data/references/{genome}/ref.fasta.{suffix}",
            genome=reference_build,
            suffix=["ann", "amb", "bwt.2bit.64", "fai", "pac", "sa", "dict"],
        ),
    output:
        bps="results/svaba/{projectid}/{sampleid}.svaba.bps.txt.gz",
        contigs="results/svaba/{projectid}/{sampleid}.svaba.contigs.bam",
        discordants="results/svaba/{projectid}/{sampleid}.svaba.discordants.txt.gz",
        log="results/svaba/{projectid}/{sampleid}.svaba.log",
        alignments="results/svaba/{projectid}/{sampleid}.svaba.alignments.txt.gz",
        vcf_indel_unfiltered="results/svaba/{projectid}/{sampleid}.svaba.unfiltered.indel.vcf",
        vcf_sv_unfiltered="results/svaba/{projectid}/{sampleid}.svaba.unfiltered.sv.vcf",
        vcf_indel_filtered="results/svaba/{projectid}/{sampleid}.svaba.indel.vcf",
        vcf_sv_filtered="results/svaba/{projectid}/{sampleid}.svaba.sv.vcf",
    benchmark:
        "results/performance_benchmarks/svaba/{projectid}/{sampleid}.svaba.log"
    params:
        outprefix="results/svaba/{projectid}/{sampleid}.svaba",
    conda:
        "../envs/svaba.yaml"
    threads: 8
    resources:
        mem_mb="16000",
        qname="large",
    shell:
        "svaba run -p {threads} -G {input.bwa_fasta} -I -L 6 -t {input.bam} -a {params.outprefix}"
