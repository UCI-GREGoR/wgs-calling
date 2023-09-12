rule compare_sv_callers:
    """
    Generate an html report summarizing the apparent relationships
    between variation in different SV callers' output vcfs
    """
    input:
        sv_source_data=lambda wildcards: expand(
            "results/reports/sv_data/{{projectid}}/{sampleid}.sv.{merge_tool}-raw.tsv",
            sampleid=list(
                set(manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"])
            ),
            merge_tool=config["behaviors"]["sv-endpoints"][wildcards.endpoint]["sv-ensemble"][
            "merge-tool"
            ],
        ),
        r_resources="workflow/scripts/compare_sv_callers.R",
    output:
        "results/reports/{projectid}/{endpoint}/sv_caller_comparison.html",
    params:
        sv_callers=lambda wildcards: config["behaviors"]["sv-endpoints"][wildcards.endpoint][
            "sv-callers"
        ],
    conda:
        "../envs/r.yaml" if not use_containers else None
    container:
        "{}/r.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["r"]["threads"]
    resources:
        mem_mb=config_resources["r"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["r"]["queue"], config_resources["queues"]
        ),
    script:
        "../scripts/compare_sv_callers.Rmd"
