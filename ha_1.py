
import functools
import math
from matplotlib import pyplot as plt
import numpy as np

# Aufgabe 3 a

def lagrange(x, nodes, k):
    """
    Calculates the k-th lagrange-polynome for a given x.

    @type x: float
    @param x: the x value
    @type nodes: list
    @param nodes: a list of supporting points
    @type k: int
    @param k: the k-th lagrange polynome to calculate

    @rtype: float
    @returns: The value of the k-th lagrange-polynome at the point x
    """

    product = 1
    xk = nodes[k]

    for i, node in enumerate(nodes):
        if i == k:
            # we skip this loop-step to make sure,
            # that the denominator of the fraction below never becomes 0
            continue
        product *= ((x - node) / (xk - node))
    
    return product


# Aufgabe 3 b

def lagrange_interpolation(x, nodes, function_values):
    """
    Calculates the lagrange interpolation for each given function value.

    @param x: the x value
    @type x: float
    @param nodes: a list of supporting points
    @type nodes: list
    @param function_values: a list of function values, that have been calculated for each supporting points.
                            note, that len(nodes) == len(function_values) must be true
    @type function_values: list

    @rtype: float
    @return the result of the formula of the lagrange-iterpolation algorithm. 
    """
    # we want to make sure our algorithm meets this condition
    assert len(nodes) == len(function_values)
    
    return functools.reduce(lambda a, b: a + (b * lagrange(x, nodes, function_values.index(b))), function_values)

    # more verbose way of writing this algorithm without a reduce function
    # sum = 0
    # for k, value in enumerate(function_values):
    #     sum += (value * lagrange(x, nodes, k))
    
    # return sum

# Aufgabe 3 c

def plot(n):
    """
    Plots graphs for given supporting points 1...n
    """
    rowspan = 2
    colspan = 3

    x = np.arange(0, math.pi, 0.02)
    for i in range(1, (n + 1)):
        # we want to arrange different subplots listed on top of each other. One plot for n = 1, one for n = 2, ...
        subplot = plt.subplot2grid(((n * 2) + 1, colspan), ((i - 1) * rowspan, 0), rowspan=rowspan, colspan=colspan)
        
        # since the calculation of the function values already need the list of supporting points, we
        # return the supporting points and their function values.
        # Therefore we make use of the tuple data structure here, because we don't want 
        # to call _getsupporting_points twice 
        # when calling the langrange_interpolation function.
        tuple_points = _get_supporting_points_and_function_values(math.sin, i)
        supporting_points = tuple_points[0]
        function_values = tuple_points[1]
        y = map(lambda xi: lagrange_interpolation(xi, supporting_points, function_values), x)

        # Possible equivalent implementation of this algorithm
        # y = []
        # tuple_points = _get_supporting_points_and_function_values(math.sin, n)
        # supporting_points = tuple_points[0]
        # function_values = tuple_points[1]
        # for xi in x:
        #     y.append(lagrange_interpolation(xi, supporting_points, function_values))

        subplot.plot(x, y)
        subplot.set_title('Plot for n = ' + str(i))
    plt.tight_layout()
    plt.show()


def _get_supporting_points_and_function_values(function, n):
    """
    Calculates the supporting points and their function values for given supporting points.

    @param function: any function that has an input paramter and that calculates something e.g. sin(x)
    @param n: the number of supporting points.
    
    @return: a tuple containing two lists. The first entry is a list of supporting points, the second entry maps a function for each supporting point.

    @see _get_supporting_points(n)
    """
    supporting_points = _get_supporting_points(n)
    return (supporting_points, map(lambda xi: function(xi), supporting_points))

    # Possible different way of the same algorithm:

    # supporting_points = _get_supporting_points(n)
    # function_values = []
    # for item in supporting_points:
    #   function_values.append(function(item))
    # return (supporting_points, function_values)


def _get_supporting_points(n):
    """
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