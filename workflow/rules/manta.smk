rule manta_configure:
    """
    Configure manta SV caller on a single bamfile.
    """
    input:
        bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        bai="results/aligned_bams/{projectid}/{sampleid}.bai",
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
        temp("{}/manta_workdir/{{projectid}}/{{sampleid}}/runWorkflow.py".format(tempDir)),
    benchmark:
        "results/performance_benchmarks/manta_configure/{projectid}/{sampleid}.tsv"
    params:
        tmpdir=expand("{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}", tempdir=tempDir),
    conda:
        "../envs/manta.yaml" if not use_containers else None
    container:
        "{}/manta.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
        tmpdir=lambda wildcards: "{}/manta_workdir/{}/{}".format(
            tempDir, wildcards.projectid, wildcards.sampleid
        ),
    shell:
        "configManta.py --config {input.manta_config} --bam {input.bam} --reference {input.fasta} "
        "--callRegions {input.calling_region} --runDir {params.tmpdir}"


rule manta_run:
    """
    After run configuration, actually run manta SV caller.
    """
    input:
        bam="results/aligned_bams/{projectid}/{sampleid}.bam",
        bai="results/aligned_bams/{projectid}/{sampleid}.bai",
        fasta="reference_data/{}/{}/ref.fasta".format(
            config["behaviors"]["aligner"], reference_build
        ),
        fai="reference_data/{}/{}/ref.fasta.fai".format(
            config["behaviors"]["aligner"], reference_build
        ),
        script=expand(
            "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/runWorkflow.py", tempdir=tempDir
        ),
    output:
        diploid_vcf=temp(
            expand(
                "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/diploidSV.vcf.gz",
                tempdir=tempDir,
            )
        ),
        diploid_tbi=temp(
            expand(
                "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/diploidSV.vcf.gz.tbi",
                tempdir=tempDir,
            )
        ),
        candidatesv_vcf=temp(
            expand(
                "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/candidateSV.vcf.gz",
                tempdir=tempDir,
            )
        ),
        candidatesv_tbi=temp(
            expand(
                "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/candidateSV.vcf.gz.tbi",
                tempdir=tempDir,
            )
        ),
        candidatesmallindels_vcf=temp(
            expand(
                "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/candidateSmallIndels.vcf.gz",
                tempdir=tempDir,
            )
        ),
        candidatesmallindels_tbi=temp(
            expand(
                "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/candidateSmallIndels.vcf.gz.tbi",
                tempdir=tempDir,
            )
        ),
    benchmark:
        "results/performance_benchmarks/manta_run/{projectid}/{sampleid}.tsv"
    params:
        tmpdir=expand("{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}", tempdir=tempDir),
    conda:
        "../envs/manta.yaml" if not use_containers else None
    container:
        "{}/manta.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["manta"]["threads"]
    resources:
        mem_mb=config_resources["manta"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["manta"]["queue"], config_resources["queues"]
        ),
        tmpdir=lambda wildcards: "{}/manta_workdir/{}/{}".format(
            tempDir, wildcards.projectid, wildcards.sampleid
        ),
    shell:
        "python2 {input.script} -j {threads}"


rule manta_sort_output:
    """
    After running manta, sort the vcf output.
    """
    input:
        vcf=expand(
            "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/diploidSV.vcf.gz",
            tempdir=tempDir,
        ),
        tbi=expand(
            "{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}/results/variants/diploidSV.vcf.gz.tbi",
            tempdir=tempDir,
        ),
    output:
        vcf="results/manta/{projectid}/{sampleid}.manta.vcf.gz",
    benchmark:
        "results/performance_benchmarks/manta_sort_output/{projectid}/{sampleid}.tsv"
    params:
        tmpdir=expand("{tempdir}/manta_workdir/{{projectid}}/{{sampleid}}", tempdir=tempDir),
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
        tmpdir=lambda wildcards: "{}/manta_workdir/{}/{}".format(
            tempDir, wildcards.projectid, wildcards.sampleid
        ),
    shell:
        "bcftools sort --temp-dir {params.tmpdir} -O z -o {output.vcf} {input.vcf}"
