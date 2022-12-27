# inputs are harmonized, duphold-applied, duphold- and pass-filtered files
# and also collapsed???  Depends on ensembler.
# run jasmine or truvari and then apply a voting schematic of some kind (or let the ensembler do it?)
# how do jasmine/truvari handle BNDs?


rule compare_sv_callers:
    """
    Generate an html report summarizing the apparent relationships
    between variation in different SV callers' output vcfs
    """
    input:
        sv_source_data=lambda wildcards: expand(
            "results/reports/sv_data/{{projectid}}/{sampleid}.sv.svdb-raw.tsv",
            sampleid=list(
                set(manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"])
            ),
        ),
        r_resources="workflow/scripts/compare_sv_callers.R",
    output:
        "results/reports/{projectid}/sv_caller_comparison.html",
    params:
        sv_callers=config["behaviors"]["sv-callers"],
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    script:
        "../scripts/compare_sv_callers.Rmd"
