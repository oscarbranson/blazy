def issubset(target, reference):
    """
    Returns True is all elements of target are present in reference.

    Both inputs must be iterable.
    """
    return set(target).issubset(set(reference))