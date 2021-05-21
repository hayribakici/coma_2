# Berechnen Sie nun für n= 1,...,500 den Fehler der Riemann-Summe, 
# und plotten Sie diesen Fehler in einer logarithmischen Skala gegen n. 
# Werten Sie dazu einmal die Funktion an den Anfangspunkten der Teilintervalle 
# aus (d.h.q= 0), ein anderes Mal an deren Mittelpunkten (d.h. q= 0.5). 
# Vergleichen Sie Ihre Ergebnisse. Was beobachten Sie?

from matplotlib import pyplot as plt
from scipy import special

import math
import numpy as np

f1 = lambda x: math.e ** (-x ** 2)

def riemann(I,f,n,q):
    """
        Berechnet die Riemann-Summe über einem bestimmten Intervall mit
        einer variablen Anzahl an Teilintervallen.

        Parameter:
        @param I:   Integrationsintervall
        @type I:    tuple
        @param f:   Funktion
        @param n:   Anzahl der Teilintervalle
        @type n:    int
        @param q:   einen Wert, der durch ξ_k=x_{k−1}+q(x_k−x_k−1) die 
                    Lage des Wertes ξ_k festlegt. Es muss gelten
                    0 ≤ q ≤ 1
        @type q:    float

        @rtype:     float
        @return     die Riemannsumme
    """
    
    xvalues = np.linspace(I[0], I[1], n+1)
    qvalues = _get_xi_k_values(xvalues, q)

    return _get_sum(f, xvalues, qvalues)

def _get_xi_k_values(list, q):
    """
        Berechnet die Werte, die durch ξ_k=x_{k−1}+q(x_k−x_k−1) die 
        Lage des Wertes ξ_k festlegen.

        Parameter:
        @param list:    Liste der x-Werte
        @type list:        list
        @param q:       Benötigt zur Berechnung von xi. Es muss gelten
                        0 ≤ q ≤ 1
        @type q:        float

        @rtype:         list
        @return         Eine Liste von xi-Werten
    """
    qv = []

    for i in range(len(list)-1):
        qv.append(list[i]+q*(list[i+1]-list[i]))
    
    return qv

def _get_sum(f, x, q):
    """
        Berechnet die Riemann-Summe.

        Parameter:
        @param f:   Funktion
        @param x:   Eine Liste der x-Werte
        @type x:    list
        @param q:   Eine Liste der xi-Werte
        @type q:    list

        @rtype:     float
        @return     Die Riemann-Summe
    """
    result = 0

    for i in range(len(q)):
        result += f(q[i])*(x[i+1]-x[i])

    return result

def _make_plot_list_with_riemann_failure(I, f, n, q):
    """
        Erzeugt zwei Listen zum plotten der Riemann-Fehler
        in einem gegebenen Intervall.

        Parameter:
        @param I:   Integrationsintervall
        @type I:    tuple
        @param f:   Funktion
        @param n:   Anzahl der Teilintervalle
        @type n:    int
        @param q:   einen Wert, der durch ξ_k=x_{k−1}+q(x_k−x_k−1) die 
                    Lage des Wertes ξ_k festlegt. Es muss gelten
                    0 ≤ q ≤ 1
        @type q:    float

        @rtype:     tupel of lists
        @return     Ein Tupel der Form (list, list)
    """
    x = list(np.arange(1,n+1,1))
    y = []
    for i in x:
        # Fehler der Riemann-Summe in Bezug zum Vergleichswert
        # 0.5*scipy.special.erf(1)*math.sqrt(math.pi) 
        y.append(riemann(I, f, i, q)-0.5*special.erf(1)*math.sqrt(math.pi))
    
    return (x, y)

def make_plot(I, f, n, q, c):

    x = _make_plot_list_with_riemann_failure(I, f, n, q)[0]
    y = _make_plot_list_with_riemann_failure(I, f, n, q)[1]

    plt.plot(x, y, c, label='q=%s' % q) 

I = (0, 1)
n = 500
# q = 0 für left-points
# q=0.5 für mid-points

make_plot(I, f1, n, 0, 'r')
make_plot(I, f1, n, 0.5, 'b')

plt.xscale("log")
plt.grid(which='both')
plt.title('Fehler der Riemann-Summe für\n left-points(q=0) und mid-points(q=0.5)')
plt.legend()
plt.xlabel('n=%i (logarithmische Skala)' %n)
plt.show()