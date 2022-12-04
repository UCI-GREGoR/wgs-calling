rule somalier_extract:
    """
    Run somalier extract on a single bam.
    """
    input:
        bam="results/markdups/{projectid}/{fileprefix}.mrkdup.sort.bam",
        bai="results/markdups/{projectid}/{fileprefix}.mrkdup.sort.bam.bai",
        fasta="reference_data/references/{}/ref.fasta".format(reference_build),
        fai="reference_data/references/{}/ref.fasta.fai".format(reference_build),
    output:
        "results/somalier/{projectid}/extract/{fileprefix}.somalier",
    benchmark:
        "results/performance_benchmarks/somalier_extract/{projectid}/{fileprefix}.tsv"
    params:
        extract_dir="results/somalier/{projectid}/extract",
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "somalier extract -d {params.extract_dir} "
        "--sites ${{CONDA_PREFIX}}/share/somalier/sites.GRCh37.vcf.gz "
        "-f {input.fasta} "
        "{input.bam}"


rule somalier_relate:
    """
    Compute relatedness metrics on preprocessed alignment data with somalier.
    """
    input:
        somalier=lambda wildcards: tc.construct_somalier_extract_targets(wildcards, manifest),
        ped="results/somalier/{projectid}/relate/somalier.ped",
    output:
        html="results/somalier/{projectid}/relate/somalier.html",
        pairs="results/somalier/{projectid}/relate/somalier.pairs.tsv",
        samples="results/somalier/{projectid}/relate/somalier.samples.tsv",
    benchmark:
        "results/performance_benchmarks/somalier_relate/{projectid}.tsv"
    params:
        outprefix="results/somalier/{projectid}/relate/somalier",
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "somalier relate --ped {input.ped} -o {params.outprefix} {input.somalier}"


rule somalier_build_pedfile:
    """
    Generate a pedfile for sex check for somalier.

    In v0.0.1 here, this pedfile will only contain placeholders.
    """
    output:
        "results/somalier/{projectid}/relate/somalier.ped",
    benchmark:
        "results/performance_benchmarks/somalier/{projectid}/relate/somalier.tsv"
    params:
        subjectids=lambda wildcards: manifest.loc[
            manifest["projectid"] == wildcards.projectid, "sampleid"
        ].to_list(),
    script:
        "../scripts/construct_somalier_pedfile.py"
