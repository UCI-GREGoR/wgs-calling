from warnings import warn


def format_reference_build(ref: str) -> str:
    """
    Given a config key for a reference genome, try to turn it
    into a formatted string describing the genome
    """
    if ref.lower() == "grch38" or ref.lower() == "hg38":
        return "GRCh38"
    if ref.lower() == "grch37" or ref.lower() == "hg19":
        return "GRCh37"
    warn(
        "The user-provided genome build configuration key '"
        + ref
        + "' does not match any of the simple expected patterns "
        "(grch38, hg38, grch37, hg19), and as such will be recapitulated "
        "as-is; this may cause problems with some tool compatibility if "
        "a tool sniffs the VCF header.",
        UserWarning,
    )
    return ref
