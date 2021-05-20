
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

def max_norm(function, interpolator, a, b, nodes):
    n = 1000
    results = []
    for i in range(n + 1):
        x = _get_supporting_point(i, a, b, n)
        results.append(abs(function(x) - interpolator(x, function, nodes)))
    return max(results)

def plot(a, b, n, interpolator_type):
    max_norms = []
    x = []
    for k in range(1, n + 1): 
        x.append(k)
        nodes = _get_supporting_points(a, b, k)
        interpolator = None
        
        if interpolator_type == InterpolationType.NEWTON:
            interpolator = newton
        else:
            interpolator = lagrange
        
        max_norms.append(max_norm(f4, interpolator, a, b, nodes))

    plt.plot(x, [f3(xi) for xi in x], color="blue")
    plt.plot(x, max_norm)
    plt.show()

def plotf():
    rowspan = 4
    colspan = 1
    x = np.arange(-5, 5, 0.02)

    for i, f in enumerate(functions):
        # we want to arrange different subplots listed on top of each other. 
        # One plot for n = 1, one for n = 2, ...
        subplot = plt.subplot(rowspan, colspan, i + 1)

        supporting_points = _get_supporting_points(-5, 5, 10)
    
        newtonY = list(map(lambda xi: newton(xi, f, supporting_points), x))
    
        subplot.plot(x, newtonY, color="blue")
        subplot.plot(x, [f(xi) for xi in x], color="red")
        subplot.set_title('Plot functions[' + str(i) + ']')
    plt.show()
        
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
        coeff_helper[i] = []
        coeff_helper[i].append(function(node))
        
        for k in range(i):
            coeff_helper[i].append((coeff_helper[i - 1][k] - coeff_helper[i][k]) / (nodes[i - k - 1] - node)) 
        
        coefficients.append(coeff_helper[i][len(coeff_helper[i]) - 1])
    
    return coefficients

def _newton_interpolation(x, coeff, nodes):
    """
    Calculates the Newton Polynome interpolation at a point x with given coefficiants and
    nodes.
    """
    result = coeff[0]
    for i in range(1, len(coeff)):
        product = coeff[i]
        for k in range(i):
            product *= (x - nodes[k])
        result += product
    return result

def newton(x, function, nodes):
    return _newton_interpolation(x, _calc_coefficients(function, nodes), nodes)

def _lagrange_k(x, nodes, k):
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

def lagrange(x, function, nodes):
    """
    Calculates the lagrange interpolation for each given function value.

    @param x: the x value
    @type x: float
    @param nodes: a list of supporting points
    @type nodes: list
    @param function: the function to apply to the nodes on

    @rtype: float
    @return the result of the formula of the lagrange-iterpolation algorithm. 
    """
    function_values = [function(node) for node in nodes]
    
    return functools.reduce(lambda a, b: a + (b * _lagrange_k(x, 
                                                    nodes,
                                                    function_values.index(b))
                                             ), function_values)

class InterpolationType(Enum):
    LAGRANGE = 1,
    NEWTON = 2

plot(-5, 5, 10, InterpolationType.NEWTON)