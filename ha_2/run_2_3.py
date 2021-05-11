import math
from matplotlib import pyplot as plt
import numpy as np


class AitkenNeville():
    """
    Class that implements the Aiken-Neville algorithm to interpolate the difference quotient at a point x.
    """
    def __init__(self, function, x, xs):
        """
        Constructor.
        @param function: the function to interpolate
        @param x: the value to interpolate
        @type x: float

        @param xs:
        @type xs: float
        """
        self._nodes = []
        self._results = {}
        self._function = function
        self._x = x
        self._xs = xs

    def difference_quotient(self, h):
        if h == 0:
            return 0
        return (self._function(self._xs + h) - self._function(self._xs)) / h

    def add_nodes(self, nodes):
        lenBefore = len(self._nodes)
        self._nodes.extend(nodes)

        self._aitken(lenBefore)

    def add_node(self, node):
        """
        Adds a node and calculates its interpolation respectively and fills its values
        into C{results}.
        @param node: the node to add. Must be the smallest number. 
        
        """
        # Checking for emptyness and out new value fulfills the requirement
        # as being the smallest item
        if self._nodes and min(self._nodes) < node:
            raise ValueError("New node must be smaller than the rest")

        lenBefore = len(self._nodes)
        
        self._nodes.append(node)
        self._aitken(lenBefore)

       
    def _aitken(self, lenBefore):
        nodes = self._nodes
        results = self._results
        x = self._x

        # We want to fill the results dictionary as with its polynome values as follows:
        # {
        #   0: [p00], 
        #   1: [p11, p01], 
        #   2: [p22, p12, p02], 
        #   3: [p33, p23, p13, p03], 
        #   4: [p44, p34, p24, p14, p04], 
        #   5: ...
        # }
        # Adding another node uses the already existing calculations.
        for k in range(lenBefore, len(self._nodes)):
            results[k] = []
            # we count backwards: from k to 0.
            for i in range(k, -1, -1):
                if i == k:
                    results[k].append(self.difference_quotient(h=nodes[k]))
                    continue

                xi = nodes[i]
                xk = nodes[k]
                frac = 1 / (xk - xi)
                # using the "p_(i+1)k - th polynome"
                p1 = (x - xi) * results[k][len(results[k]) - 1]
                # using the "p_i(k-1) - th polynome"
                p2 = (x - xk) * results[k - 1][(len(results[k - 1]) - 1) - i]

                # putting everythin together
                results[k].append(frac * (p1 - p2))


    def get_polynome(self, n): 
        """
        @param n: the n-th polynome to retrieve.
        @return ... TODO
        """
        polynome = len(self._results[n]) - 1
        return self._results[n][polynome]

def f_3(x): 
    gamma = 10000
    return (1 / gamma) * math.atan(gamma * x)

def test(n):
    x_ = 0
    H = math.pi / 2

    
    nodes_1 = _get_supporting_points(n, lambda k: 2**(-k) * H)
    nodes_2 = _get_supporting_points(n, lambda k: H / (k + 1)) # have length of n



    # 1. f'(x)
    a = AitkenNeville(lambda x: 1 + math.sin(x), 0, 0)
    a.add_nodes(nodes_1)
    pn = a.get_polynome(n)

    error_wo_extra = abs(math.cos(x_) - a.difference_quotient(nodes_1[len(nodes_1) - 1]))
    error_w_extra = abs(math.cos(x_) - pn)
    print(error_wo_extra)
    print(error_w_extra)


    # n = 1, ..., 40
    # hk = 2^-k * pi / 2, k = 0, ..., n
    # for i in range(n):    

        # a = AitkenNeville(lambda x: 1 + math.sin(x), 0, 0)
        
        # b = AitkenNeville(lambda x: math.sqrt(x**3), 0, 0)

        # c = AitkenNeville(f_3, 0, 0)
    
def _get_supporting_points(n, function):
    mylist = []
    for k in range(n + 1): # 0, ..., n
        mylist.append(function(k))

    return mylist

test(n=40)