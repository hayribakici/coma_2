import math
from matplotlib import pyplot as plt
import numpy as np

def difference_quotient(h, function, x):
    if h == 0:
        return 0
    return (function(x + h) - function(x)) / h


def aitken_neville(x, function, nodes):
    n = len(nodes) - 1
    return _aitken_neville(0, n, function, x, nodes)

def _aitken_neville(i, k, function, x, nodes):
    xi = nodes[i]
    xk = nodes[k]
    print("i: %s, k: %s, xi: %s, xk: %s", i, k, xi, xk)
    if i == k:
        f = difference_quotient(xi, function, 1)
        print(f)
        return f
    
    return (1 / (xk - xi)) * (((x - xi) * _aitken_neville(i + 1, k, function, x, nodes) - 
                (x - xk) * _aitken_neville(i, k - 1, function, x, nodes)))
    

def test():
    nodes = [0.2, 0.4] #np.arange(0, 1, 0.2)
    print(nodes)
    print(aitken_neville(0, math.log, nodes))
    nodes.append(0.6)
    print(aitken_neville(0, math.log, nodes))

test()