rule somalier_extract:
    """
    Run somalier extract on a single bam.
    """
    input:
        bam="{pathprefix}/markdups/{projectid}/{fileprefix}.mrkdup.sort.bam",
        bai="{pathprefix}/markdups/{projectid}/{fileprefix}.mrkdup.sort.bam.bai",
        fasta=config["references"]["grch37"]["fasta"],
    output:
        "{pathprefix}/somalier/{projectid}/extract/{fileprefix}.somalier",
    params:
        extract_dir="{pathprefix}/somalier/{projectid}/extract",
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        h_vmem="2000",
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
        somalier=lambda wildcards: tc.construct_somalier_extract_targets(manifest),
        ped="{pathprefix}/somalier/{projectid}/relate/somalier.ped",
    output:
        html="{pathprefix}/somalier/{projectid}/relate/somalier.html",
        pairs="{pathprefix}/somalier/{projectid}/relate/somalier.pairs.tsv",
        samples="{pathprefix}/somalier/{projectid}/relate/somalier.samples.tsv",
    params:
        outprefix="{pathprefix}/somalier/{projectid}/relate/somalier",
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        h_vmem="4000",
        qname="small",
    shell:
        "somalier relate --ped {input.ped} -o {params.outprefix} {input.somalier}"


rule somalier_build_pedfile:
    """
    Generate a pedfile for sex check for somalier.

    In v0.0.1 here, this pedfile will only contain placeholders.
    """
    output:
        "{pathprefix}/somalier/{projectid}/relate/somalier.ped",
    params:
        subjectids=lambda wildcards: manifest.loc[
            manifest["projectid"] == wildcards.projectid, "sampleid"
        ].to_list(),
    script:
        "../scripts/construct_somalier_pedfile.py"
