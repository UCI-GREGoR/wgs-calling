from warnings import warn


def format_reference_build(ref: str) -> str:
    """
    Given a config key for a reference genome, try to turn it
    into a string that Moon will acknowledge

    So here's the thing: Moon is a very strange tool that implicitly
    enforces a bunch of restrictions on otherwise valid vcf and then
    doesn't bother to tell you about it. One of those restrictions is
    for its support of GRCh38: it tries to determine reference genome
    with the following logic:

    - if it finds a vcf header entry that contains anywhere in the line
    '##reference=', it takes the rest of the line after those characters
    and searches exactly for 'hg19' or 'GRCh37'; and then if that fails
    it searches for 'hg38' or 'GRCh38'
    - if the above fails for any reason, it then sniffs the entire header
    first for 'GRCh38' or 'hg38'
    - if the above fails, it returns 37

    Note that the '##reference' search is not anchored, and the 'GRCh' searches
    are case-sensitive
    """
    if ref.lower() == "grch38" or ref.lower() == "hg38":
        return "GRCh38"
    if ref.lower() == "grch37" or ref.lower() == "hg19":
        return "GRCh37"
    warn(
        "The user-provided genome build configuration key '"
        + ref
        + "' does not match any of the simple patterns expected by Moon; "
        "as such, this will cause mismatch problems if you try to load "
        "a variant VCF into Moon without first modifying the header to "
        "contain 'GRCh38' or 'hg38' somewhere in the header. If your "
        "intended build is GRCh37/hg19, Moon should default to this anyway."
    )
    return ref
