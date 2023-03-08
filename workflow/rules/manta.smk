rule manta_configure:
    """
    Configure manta SV caller on a single bamfile.
    """
    input:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        manta_config=config["parameters"]["manta"]["config-ini"],
        calling_region="reference_data/manta/{}/ref.calling.range.bed.gz".format(reference_build),
        calling_region_tbi="reference_data/manta/{}/ref.calling.range.bed.gz.tbi".format(
            reference_build
        ),
    output:
        temp("temp/manta_workdir/{projectid}/{sampleid}/runWorkflow.py"),
    benchmark:
        "results/performance_benchmarks/manta_configure/{projectid}/{sampleid}.tsv"
    params:
        tmpdir="temp/manta_workdir/{projectid}/{sampleid}",
    conda:
        "../envs/manta.yaml"
    container:
        "apptainer_images/manta.sif"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
        tmpdir=lambda wildcards: "temp/manta_workdir/{}/{}".format(
            wildcards.projectid, wildcards.sampleid
        ),
    shell:
        "configManta.py --config {input.manta_config} --bam {input.bam} --reference {input.fasta} "
        "--callRegions {input.calling_region} --runDir {params.tmpdir}"


rule manta_run:
    """
    After run configuration, actually run manta SV caller.
    """
    input:
        bam="results/bqsr/{projectid}/{sampleid}.bam",
        bai="results/bqsr/{projectid}/{sampleid}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        script="temp/manta_workdir/{projectid}/{sampleid}/runWorkflow.py",
    output:
        diploid_vcf=temp(
            "temp/manta_workdir/{projectid}/{sampleid}/results/variants/diploidSV.vcf.gz"
        ),
        diploid_tbi=temp(
            "temp/manta_workdir/{projectid}/{sampleid}/results/variants/diploidSV.vcf.gz.tbi"
        ),
        candidatesv_vcf=temp(
            "temp/manta_workdir/{projectid}/{sampleid}/results/variants/candidateSV.vcf.gz"
        ),
        candidatesv_tbi=temp(
            "temp/manta_workdir/{projectid}/{sampleid}/results/variants/candidateSV.vcf.gz.tbi"
        ),
        candidatesmallindels_vcf=temp(
            "temp/manta_workdir/{projectid}/{sampleid}/results/variants/candidateSmallIndels.vcf.gz"
        ),
        candidatesmallindels_tbi=temp(
            "temp/manta_workdir/{projectid}/{sampleid}/results/variants/candidateSmallIndels.vcf.gz.tbi"
        ),
    benchmark:
        "results/performance_benchmarks/manta_run/{projectid}/{sampleid}.tsv"
    params:
        tmpdir="temp/manta_workdir/{projectid}/{sampleid}",
    conda:
        "../envs/manta.yaml"
    container:
        "apptainer_images/manta.sif"
    threads: 4
    resources:
        mem_mb="24000",
        qname="small",
        tmpdir=lambda wildcards: "temp/manta_workdir/{}/{}".format(
            wildcards.projectid, wildcards.sampleid
        ),
    shell:
        "python2 {input.script} -j {threads}"


rule manta_sort_output:
    """
    After running manta, sort the vcf output.
    """
    input:
        vcf="temp/manta_workdir/{projectid}/{sampleid}/results/variants/diploidSV.vcf.gz",
        tbi="temp/manta_workdir/{projectid}/{sampleid}/results/variants/diploidSV.vcf.gz.tbi",
    output:
        vcf="results/manta/{projectid}/{sampleid}.manta.vcf.gz",
    benchmark:
        "results/performance_benchmarks/manta_sort_output/{projectid}/{sampleid}.tsv"
    params:
        tmpdir="temp/manta_workdir/{projectid}/{sampleid}",
    conda:
        "../envs/bcftools.yaml"
    container:
        "apptainer_images/bcftools.sif"
    threads: 4
    resources:
        mem_mb="16000",
        qname="small",
        tmpdir=lambda wildcards: "temp/manta_workdir/{}/{}".format(
            wildcards.projectid, wildcards.sampleid
        ),
    shell:
        "bcftools sort --temp-dir {params.tmpdir} -O z -o {output.vcf} {input.vcf}"
