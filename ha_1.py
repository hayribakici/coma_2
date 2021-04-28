

# Aufgabe 3 a

def lagrange(x, nodes, k):
    """
    Calculates the lagrange-polynome.

    @type x: float
    @param x: the x value
    @param nodes: a list of supporting points
    @type k: int
    @param k: the k-th lagrange polynome to calculate

    @rtype: float
    @returns: The value of the lagrange-polynome at the point x
    """

    product = 1
    for i, node in enumerate(nodes):
        xk = 1 if i == k else 0
        product *= ((x - node) / (xk - node))
    
    return product


# Aufgabe 3 b

def lagrange_interpolation(x, nodes, function_values):
    # TODO implement me
    pass

# Aufgabe 3 c

def plot():
    # TODO implement me
    pass
