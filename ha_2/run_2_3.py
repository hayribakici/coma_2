import math
from matplotlib import pyplot as plt
import numpy as np

def difference_quotient(h, function, x):
    if h == 0:
        return 0
    return (function(x + h) - function(x)) / h


def aitken_neville(x, x_, function, nodes):
    n = len(nodes) - 1
    return _aitken_neville(0, n, x, x_, function, nodes)

def _aitken_neville(i, k, x, x_, function, nodes):
    xi = nodes[i]
    xk = nodes[k]
    print("p_", i, "_", k, " = 1 / (", xk, "-", xi, ") * ((", x, "-", xi, ") * p_", i + 1, "_", k, " - (", x, " - ", xk, ") * p_", i, "_", k - 1, ")")  
    
    if i == k:
        dx = difference_quotient(xi, function, x_)
        print(i, " == ", k, " d/dx = ", dx)
        return dx   
    
    return (1 / (xk - xi)) * (((x - xi) * _aitken_neville(i + 1, k, x, x_, function, nodes)) - 
                ((x - xk) * _aitken_neville(i, k - 1, x, x_, function, nodes)))
    

def test():
    nodes = [0.2, 0.4, 0.6] #np.arange(0, 1, 0.2)
    print(aitken_neville(0, 1, math.log, nodes))
    # TODO what to when we add another node?
    # nodes.append(0.6)
    # print(_aitken_neville(len(nodes) - 2, len(nodes) - 1, math.log, 0, nodes))

test()