
import functools
import itertools
import math
from matplotlib import pyplot as plt
import numpy as np

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

    # return functools.reduce(lambda a, xi: a * (x - xi) / (nodes[k] - xi), nodes)

    # Possible different way of writing this algorithm.
    # TODO this throws an exception when nodes[k] - node = 0
    # looks like nodes[k] (or the x_k to be exact) is a different number
    product = 1
    for node in nodes:
        product *= ((x - node) / (nodes[k] - node))
    
    return product


# Aufgabe 3 b

def lagrange_interpolation(x, nodes, function_values):
    """
    Calculates the lagrange interpolation.
    """
    # we want to make sure our algorithm meets this condition
    assert len(nodes) == len(function_values)
    
    # return functools.reduce(lambda a, b: a + b * lagrange(x, nodes, function_values.index(b)), function_values)

    sum = 0
    for k, value in enumerate(function_values):
        sum += (value * lagrange(x, nodes, k))
    
    return sum

# Aufgabe 3 c

def plot(n):
    """
    Plots in a grid of n rows
    """


    rowspan = 2
    colspan = 3
    x = np.arange(0, math.pi, 0.02)
    y = map(lambda item: lagrange_interpolation(item, _get_supporting_points(n), _get_function_values(math.sin, n)), x)
    plt.plot(x, y)

    # TODO uncomment when above code workds
    # for i in range(1, (n + 1)):
    #     subplot = plt.subplot2grid(((n * 2) + 1, colspan), ((i - 1) * rowspan, 0), rowspan=rowspan, colspan=colspan)
    #     y = map(lambda item: lagrange_interpolation(item, _get_supporting_points(i), _get_function_values(math.sin, i)), x)
    #     subplot.plot(x, y)
    #     subplot.set_title('for n = ' + str(i))
    # plt.tight_layout()
    plt.show()


def _get_function_values(function, n):
    """
    Calculates the function values for given supporting points.

    @param function: any function that has an input paramter and that calculates something e.g. sin(x)
    @param n: the number of supporting points.
    
    @return: a list of function values calculated from the supporting points

    @see _get_supporting_points(n)
    """
    return map(lambda item: function(item), _get_supporting_points(n))

    # Possible different way of the same algorithm:

    # supporting_points = _get_supporting_points(n)
    # function_values = []
    # for item in supporting_points:
    #   function_values.append(function(item))
    # return function_values


def _get_supporting_points(n):
    """
    @type k: int
    @type n: int

    @rtype: list
    @return: a list of supporting points with k * pi / n with k = 0...n
    """
    # since k has the range from 0 ... n, the range to n + 1 is necessary
    # otherwise it would go from 0 ... n - 1
    mylist = [k for k in range(n + 1)] # equals [0, 1, 2, ..., n]

    return map(lambda item: (item * math.pi) / n, mylist)

    # Possible different way of the same algorithm:

    # new_list = []
    # for item in list:
    #   new_list.append((item * pi) / n)
    # return new_list


plot(n=4)