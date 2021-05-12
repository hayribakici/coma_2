import math
from matplotlib import pyplot as plt
import numpy as np
import inspect
import sympy as sym

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
        """
        Calculates the difference quotient for a given h.
        """
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

        # We want to fill the results dictionary with its polynome values as follows:
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

class Coma2Function():
    """
    Helper class.
    """

    def __init__(self, pretty, f, ddxpretty, ddxf):
        self.pretty_function = pretty
        self.function = f
        self.ddx_pretty = ddxpretty
        self.ddxfunction = ddxf

class Test():
    def __init__(self, n, xs):
        self._n = n
        self._xs = xs
        self._H = sym.pi / 2
        self.k = sym.Symbol('k')
        self._functions =  [Coma2Function("1 + sin(x)", lambda xs: 1 + math.sin(xs), "cos(x)", lambda xs: math.cos(xs)),
                            Coma2Function("sqrt(x^3)", lambda xs: math.sqrt(xs**3), "3 * sqrt(x) / 2", lambda xs: 1.5 * math.sqrt(xs)),
                            Coma2Function("1/10000 * atan(10000 * x)", self.f_3, "1/(100000000 * x^2 + 1)", lambda xs: 1/(100000000 * xs**2 + 1))]
        self._node_functions = [2 ** (-self.k) * self._H, self._H / (self.k + 1)]
        

    def f_3(self, x): 
        gamma = 10000
        return (1 / gamma) * math.atan(gamma * x)

    def calculate_errors(self):
        """
        Calculates the error with and without extrapolation for each n.
        """
        n = self._n
        x_ = self._xs
        functions = self._functions
        

        for i in range(1, n + 1):
            print("======================= n =", i, "=======================")
            print("")

            for f in functions:

                for node_function in self._node_functions:
                    print("Current nodes:", node_function)
                    nodes = self.__get_supporting_points(i, sym.lambdify(self.k, node_function))
                    hn = self.__hn(nodes)

                    print("Calculating errors for function ", f.pretty_function)
                    neville = AitkenNeville(f.function, 0, self._xs)
                    neville.add_nodes(nodes)
                        
                    derived = f.ddxfunction

                    dq = neville.difference_quotient(hn)
                    print("|", f.ddx_pretty, " - D(", hn, ")| =", self.__get_error(derived, x_, dq))
                    print("")
                    
                    pn = neville.get_polynome(i)
                    print("|", f.ddx_pretty, " - p_", i, "| =", self.__get_error(derived, x_, pn))
                    print("")

        print("Hallo Lenni!")
        print("")
        print("Leider war uns die Aufgabenstellung nicht klar genug, sodass wir dir")
        print("noch zeitgerecht schreiben konnten.")
        print("")
        print("Uns war nicht klar, ob wir die Stützstellen h_k mit wachsendem n (von 1, ..., 40)")
        print("behandeln sollen (das würden dann insgesamt 40 plot machen!!) oder einfach")
        print("mit jedem n = n + 1 eine neue Stützselle berechnet und den Algorithmus gefüttert")
        print("werden muss (das Zweite würde eher Sinn machen). Sollte dies sein, so stellten wir")
        print("uns die weitere Frage, wie das überhaupt geplottet werden soll?")
        print("Besser gesagt, wie sollen die Fehler mit/ohne Extrapolation dargestellt werden?")
        print("Sorry, dass wir nur eine halbgare Lösung abgegeben haben, ich hoffe wir bekommen")
        print("trotzdem eine faire Punktzahl.")
        print("")
        print("Viele Grüße")
        print("Monia, Jan und Hayri")


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

def main():
    test = Test(n=40, xs=0)
    test.calculate_errors()

if __name__ == "__main__":
    main()