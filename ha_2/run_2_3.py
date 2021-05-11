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
        """
        Adds a collection of nodes to calculate its interpolation and fills its values
        into C{results}.
        @param nodes: the nodes to add.

        @see get_polynome(n) to get the n-th value
        """
        lenBefore = len(self._nodes)
        self._nodes.extend(nodes)

        self.__aitken(lenBefore)

    def add_node(self, node):
        """
        Adds a node and calculates its interpolation respectively and fills its values
        into C{results}.
        @param node: the node to add. Must be the smallest number. 
        
        @see get_polynome(n) to get the n-th value
        """
        # Checking for emptyness and out new value fulfills the requirement
        # as being the smallest item
        if self._nodes and min(self._nodes) < node:
            raise ValueError("New node must be smaller than the rest")

        lenBefore = len(self._nodes)
        
        self._nodes.append(node)
        self.__aitken(lenBefore)

       
    def __aitken(self, lenBefore):
        """
        Calculates the aitken-neville algorithm.
        """
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
        # Adding another node(s) uses the already existing calculations.
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

                # putting everything together
                results[k].append(frac * (p1 - p2))


    def get_polynome(self, n): 
        """
        @param n: the n-th polynome to retrieve.
        @return the n-th polynome that was calculated with C{add_node} or C{add_nodes} 
            or C{0} when the polynome was not found.
        """
        if not n in self._results:
            return 0
        polynome = len(self._results[n]) - 1
        return self._results[n][polynome]


class Test():

    def __init__(self, n):
        self._n = n
        self._x = 0
        self._H = math.pi / 2
        self._functions = [lambda x: 1 + math.sin(x), lambda x: math.sqrt(x**3), self.f_3]
        self._derived_functions = [lambda x: math.cos(x), lambda x: 1.5 * math.sqrt(x), lambda x: 1/(100000000* x**2 + 1)]

    def f_3(self, x): 
        gamma = 10000
        return (1 / gamma) * math.atan(gamma * x)

    def calculate(self):
        n = self._n
        x_ = self._x
        H = self._H
        

        for i in range(1, n + 1):
            nodes = self.__get_supporting_points(i, )
            self.__calculate(i, )

    def __calculate(self, i, nodes):
        functions = self._functions
        derived_functions = self._derived_functions

        hn = self.__hn(nodes)
        for f in functions:
            neville = AitkenNeville(f, 0, 0)
            neville.add_nodes(nodes)
                
            dx = neville.difference_quotient(hn)



            pn = neville.get_polynome(i)


    def __get_supporting_points(self, n, function): 
        """
        @param n: the number of nodes to return
        @param function: a function to apply for each k = 0, ..., n

        @rtype: list
        @return C{function} applied for each k = 0, ..., n
        """
        return [function(k) for k in range(n + 1)]

    def __hn(self, nodes): 
        """
        @return the last element of C{nodes}
        @rtype: float
        """
        return nodes[len(nodes) - 1]

    def __get_error(self, derived_function, x, value):
        """
        Calculates the error: |f'(x) - value|
        @param derived_function: the derived function
        @param x: the position for the derived function to calculate
        @type x: float
        @param value: the value to calculate the error.

        @return the error between C{derived_function} and C{value}
        @rtype: float
        """
        return abs(derived_function(x) - value)

    # 1. f'(x)
    a = AitkenNeville(lambda x: 1 + math.sin(x), 0, 0)
    a.add_nodes(nodes_1)
    pn = a.get_polynome(n)

    error_wo_extra = abs(math.cos(x_) - a.difference_quotient(nodes_1[len(nodes_1) - 1]))
    error_w_extra = abs(math.cos(x_) - pn)
    print(error_wo_extra)
    print(error_w_extra)
    

test(n=40)