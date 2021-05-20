# riemann(I,f,n,q)
# Vektor I: Integrationsintervall
# f: Funktion
# n: Anzahl der Teilintervalle
# 0 ≤ q ≤ 1 einen Wert, der durch ξ_k=x_{k−1}+q(x_k−x_k−1) die Lage des Wertes ξ_k festlegt

# Berechnen Sie nun für n= 1,...,500 den Fehler der Riemann-Summe, 
# und plotten Sie diesen Fehler in einer logarithmischen Skala gegen n. 
# Werten Sie dazu einmal die Funktion an den Anfangspunkten der Teilintervalle 
# aus (d.h.q= 0), ein anderes Malan deren Mittelpunkten (d.h. q= 0.5). 
# Vergleichen Sie Ihre Ergebnisse. Was beobachten Sie?

from matplotlib import pyplot as plt
from scipy import special

import math
import numpy as np



f1 = lambda x: math.e ** (-x ** 2)

def riemann(I,f,n,q):
    
    xvalues = np.linspace(I[0], I[1], n+1)

    qvalues = []

    for i in range(len(xvalues)-1):
        qvalues.append(xvalues[i]+q*(xvalues[i+1]-xvalues[i]))

    result = 0

    for i in range(len(qvalues)):
        result += f(qvalues[i])*(xvalues[i+1]-xvalues[i])

    return result

def plot(I, n, q):
    x = []
    y0 = []
    y1 = []
    for i in range(1, n):
        x.append(i)
        y0.append(riemann(I, f1, i, q)-0.5*special.erf(1)*math.sqrt(math.pi))
        y1.append(riemann(I, f1, i, 0.5)-0.5*special.erf(1)*math.sqrt(math.pi))
    
    
    plt.plot(x, y0, 'r', label='q=0')
    plt.plot(x, y1, 'b', label='q=0.5')
    plt.title('Fehler der Riemann-Summe')
    plt.legend()
    plt.xlabel('n')
    plt.show()

I = (0, 1)
n = 500
q = 0   # q=0 für left-points, q=0.5 für mid-points

plot(I, n, q)