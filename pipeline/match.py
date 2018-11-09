def match_by_fbtr(fbtr_list1, fbtr_list2):
    """Match FlyBase transcript ID's between 2 lists.

    Parameters
    ----------
    fbtr_list1 : list
    fbtr_list2 : list

    Returns
    -------
    list
        FlyBase transcript ID's found in both lists
    """
    fbtr_matches = set()
    for k in fbtr_list1:
        if k in fbtr_list2:
            fbtr_matches.add(k)
    return list(fbtr_matches)

def reciprocal_match_by_fbtr(fbtr_list1, fbtr_list2):
    """Reciprocally match FlyBase transcript ID's between 2 lists.

    Parameters
    ----------
    fbtr_list1 : list
    fbtr_list2 : list

    Returns
    -------
    list
        FlyBase transcript ID's found in both lists in both
        forward and reverse matching.
    """
    forward_match = set(match_by_fbtr(fbtr_list1, fbtr_list2))
    reverse_match = set(match_by_fbtr(fbtr_list2, fbtr_list1))
    return forward_match.intersection(reverse_match)
