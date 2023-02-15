import re


def get_svtype(chrom: str, alt: str) -> str:
    """
    Derived from https://github.com/walaj/svaba/issues/4
    """
    pattern = {"++": "INV", "-+": "INS", "+-": "DEL", "--": "INV"}
    frst = alt.split(":")[0]
    try:
        chrom2 = re.search(r"\d+", frst).group()
    except Exception:
        chrom2 = re.search(r"\d+", frst[-1])

    if frst[0].isalpha():
        res = "+"
        if frst[1] == "]":
            res += "+"
        else:
            res += "-"
    else:
        res = "-"
        if frst[0] == "]":
            res += "+"
        else:
            res += "-"

    if chrom == chrom2:
        return pattern[res]
    return "BND"


def reclassify_svs(variant: str) -> str:
    """
    Parse a bedpe line, reclassify it if needed,
    and return its output
    """
    parsed_variant = variant.split("\t")
    svtype = parsed_variant[10]
    if svtype == "BND":
        updated_svtype = get_svtype(parsed_variant[0], parsed_variant[14])
        parsed_variant[10] = updated_svtype
        parsed_variant[14] = "<{}>".format(updated_svtype)
        parsed_variant[15] = "."
        parsed_variant[16] = "."
        parsed_variant[17] = "."
        updated_info = "SVTYPE={};POS={};SVLEN={};END={}".format(
            updated_svtype,
            parsed_variant[1],
            int(parsed_variant[4]) - int(parsed_variant[1]),
            parsed_variant[4],
        )
        parsed_variant[18] = re.sub("SVTYPE=[^;]+;POS=[^;]+;", updated_info, parsed_variant[18])
        parsed_variant[19] = "."

    return "\t".join(parsed_variant)


def run_reclassify_svs(infn: str, outfn: str) -> None:
    """
    Parse an input file, and emit its reclassified SV status
    along with any necessary annotation adjustments.
    """
    variants = []
    with open(infn, "r") as f:
        variants = f.readlines()
    with open(outfn, "w") as f:
        for variant in variants:
            if variant.startswith("#"):
                f.writelines([variant])
            else:
                f.writelines([reclassify_svs(variant)])


run_reclassify_svs(
    snakemake.input[0],  # noqa: F821
    snakemake.output[0],  # noqa: F821
)
