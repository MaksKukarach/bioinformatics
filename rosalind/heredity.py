def num_pairs(n):
    return int(n * (n - 1) / 2)


def prob_dominant(dom: int, het: int, rec: int):
    total = num_pairs(dom + het + rec)
    p_both_dom = num_pairs(dom) / total
    p_both_het = (3 / 4) * (num_pairs(het) / total)
    p_dom_rec = (dom * rec) / total
    p_dom_het = (dom * het) / total
    p_het_rec = (1 / 2) * ((het * rec) / total)
    result = p_both_dom + p_both_het + p_dom_rec + p_dom_het + p_het_rec
    return round(result, 5)