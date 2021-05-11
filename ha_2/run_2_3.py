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
        self.results = {}
        self._function = function
        self._x = x
        self._xs = xs

    def __difference_quotient(self, h):
        if h == 0:
            return 0
        return (self._function(self._xs + h) - self._function(self._xs)) / h

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

        nodes = self._nodes
        lenBefore = len(self._nodes)
        
        nodes.append(node)

        results = self.results
        x = self._x

        # [p_00, p_11, p_01, p_22, p_12, p_02,     p_33, p_23, p_13, p_03]
        # p_12 = .... P_(1+1),2 - p_1,(2-1)
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
                    results[k].append(self.__difference_quotient(h=nodes[k]))
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

def test():
    longbottom = AitkenNeville(math.log, 0, 1)
    nodes = np.arange(1, 0, -0.2)
    for node in nodes:
        longbottom.add_node(node)
        print(longbottom.results)

test()