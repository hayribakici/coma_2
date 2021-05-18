
from enum import Enum
from matplotlib import pyplot as plt

import math
import numpy as np
import functools

f1 = lambda x: math.e ** x
f2 = lambda x: (x ** 2) / (1 + x ** 2)
f3 = lambda x: np.sign(x) * x ** 2
f4 = lambda x: np.sign(x) * x ** 3
functions = [f1, f2, f3, f4]

def _get_supporting_points(a, b, n):
    """
    Returns the supporting points for an interval [a, b] given n.
    @rtype list
    @return a list of supporting points with the formula: x_i = a + (i * (b - a) / n), i = 0,..., n
    """
    return [_get_supporting_point(i, a, b, n) for i in range(n + 1)]

def _get_supporting_point(i, a, b, n):
    """
    Returns a supporting point for a given i.
    """
    if i > n:
        raise ValueError("i must be < n")
    return a + ((i * (b - a)) / n)

def max_norm(values):
    return max(values)

def _get_function_values(function, nodes):
    return list(map(lambda xi: function(xi), nodes))

def plot(a, b, n, interpolation_type):
    x = np.arange(a, b, 0.02) 
    
    l = []
    
    nodes = _get_supporting_points(a, b, 1)
        
    f1(node) - newton(.., f1, nodes)
    
        
# Interpolation functions

def _calc_coefficients(function, nodes):

    # We want to fill the results dictionary with its polynome values as follows:
    # {
    #   0: y0 = [D00]
    #   1: y1 = [D01, D10]
    #   2: y2 = [D02, D11, D20]
    #   3: y3 = [D03, D12, D21 D30]
    #   4: y4 = [D04, D13, D22, D31, D40]
    #   5: ...
    # }

    coeff_helper = {}
    # holds the last element of each list: [D00, D10, D20, D30, D40, ...]
    coefficients = []
    for i, node in enumerate(nodes):
        coeff_helper[i].append(function(node))
        
        for k in range(i):
            coeff_helper[i].append((coeff_helper[i - 1][k] - coeff_helper[i][k]) / (nodes[i - k - 1] - node)) 
        
        coefficients.append(coeff_helper[i][len(coeff_helper[i]) - 1])
    
    return coefficients

def newton_interpolation(x, coeff, nodes):
    """
    Calculates the Newton Polynome interpolation at a point x with given coefficiants and
    nodes.
    """
    result = coeff[0]
    for i in range(1, len(coeff)):
        result += coeff[i] * functools.reduce(lambda a, b: (x - a) * (x - b), nodes[slice(i - 1)])

    return result

def newton(x, function, nodes):
    return newton_interpolation(x, _calc_coefficients(function, nodes), nodes)

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

def lagrange_interpolation(x, nodes, function_values):
    """
    Calculates the lagrange interpolation for each given function value.

    @param x: the x value
    @type x: float
    @param nodes: a list of supporting points
    @type nodes: list
    @param function_values: a list of function values, that have been 
                            calculated for each supporting points.
                            Note, that len(nodes) == len(function_values)
                            must be true
    @type function_values: list

    @rtype: float
    @return the result of the formula of the lagrange-iterpolation algorithm. 
    """
    # we want to make sure our algorithm meets this condition
    assert len(nodes) == len(function_values)
    
    return functools.reduce(lambda a, b: a + (b * lagrange(x, 
                                                    nodes,
                                                    function_values.index(b))
                                             ), function_values)

class InterpolationType(Enum):
    LAGRANGE = 1,
    NEWTON = 2


plot(-5, 5, 1, InterpolationType.NEWTON)