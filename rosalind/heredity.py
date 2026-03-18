def num_pairs(n):
    return int(n * (n - 1) / 2)


def prob_dominant(dom: int, het: int, rec: int):
    """Given a sample with counts of homozygous dominant, heterozygous and homozygous recessive
    organisms, return a probability that a random pair of organisms will produce and offspring
    with a dominant allele."""
    total = num_pairs(dom + het + rec)
    p_both_dom = num_pairs(dom) / total
    p_both_het = (3 / 4) * (num_pairs(het) / total)
    p_dom_rec = (dom * rec) / total
    p_dom_het = (dom * het) / total
    p_het_rec = (1 / 2) * ((het * rec) / total)
    # We use the law of total probability
    result = p_both_dom + p_both_het + p_dom_rec + p_dom_het + p_het_rec
    return round(result, 5)