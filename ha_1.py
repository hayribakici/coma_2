
import functools

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

    return functools.reduce(lambda a, xi: a * (x - xi) / (nodes[k] - xi), nodes)

    # Possible different way of writing this algorithm.
    # product = 1
    # for node in nodes:
    #     product *= ((x - node) / (nodes[k] - node))
    
    # return product


# Aufgabe 3 b

def lagrange_interpolation(x, nodes, function_values):
    """
    Calculates the lagrange-polynome ?
    """
    # we want to make sure our algorithm meets this condition
    assert len(nodes) == len(function_values)
    # TODO use indexed higher order function if possible
    functools.
    return functools.reduce(lambda a, b: a + b*lagrange(x, nodes,), function_values)


# Aufgabe 3 c

def plot():
    # TODO implement me
    pass
