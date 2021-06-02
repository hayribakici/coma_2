from matplotlib import pyplot as plt
from scipy import special

import math
import numpy as np

def summed_newton_cotes(a, b, f, k, n):
    # Seite zum Vergleich: https://de.planetcalc.com/4324/
    """
        Berechnet die summierten Newton–Cˆotes-Formel

        Parameter:
        @param a:   Startpunkt des Intervalls (a,b)
        @type a:    int
        @param b:   Endpunkt des Intervalls (a,b)
        @type b:    int
        @param f:   die zu integrierende Funktion
        @param k:   k∈{1,2,6}, gibt an, welche k-te 
                    Newton–Cˆotes-Formel angewendet wird
                    Es gilt:
                    k=1: Trapezregel    1/2  1/2
                    k=2: Simpson-Regel  1/6  4/6  1/6
                    k=6: Weddle-Regel   41/840  216/840  27/840 
                                        272/840  27/840  216/840  41/840
        @type k:    int
        @param n:   Anzahl der Intervalle, auf denen 
                    die k-te Newton–Cˆotes-Formel 
                    (mit k∈{1,2,6}) angewendet wird
        @type n:    int

        @rtype:     tupel
        @return     Ein Tupel der Form (S, A), wobei S 
                    der Wert der Quadraturformel ist und A 
                    die Anzahl der durchgeführten f-Auswertungen.
    """
    
    I = np.linspace(a, b, n+1)

    # Trapezregel
    if k==1:
        return _get_value_of_trapez_regel(I, f)
    
    # Simpson-Regel
    if k==2:
        return _get_value_of_simpson_regel(I, f)
    
    # Weddle-Regel
    if k==6:
        return _get_value_of_weddle_regel(I, f)

def _get_value_of_trapez_regel(I, f):
    S = 0
    A = 0

    for i in range(0, len(I)-1):
        a = I[i]
        b = I[i+1]
        S += (b - a) * ((f(a) + f(b)) / 2)
        A += 2

    return (S,A)

def _get_value_of_simpson_regel(I, f):
    S = 0
    A = 0

    for i in range(0, len(I) - 1):
        a = I[i]
        b = I[i + 1]
        S += ((b - a) / 6) * (f(a) + 4 * f((b + a) / 2) + f(b))
        A += 3

    return (S,A)

def _get_value_of_weddle_regel(I, f):
    S = 0
    A = 0

    for i in range(0, len(I)-1):
        a = I[i]
        b = I[i+1]
        h = (b-a)/6
        S += (((b - a)  /840) * (41 * f(a + (0 * h)) 
            + 216 * f(a + (1 * h)) + 27 * f(a + (2 * h)) 
            + 272 * f(a + (3 * h)) + 27 * f(a + (4 * h)) 
            + 216 * f(a +(5 * h)) + 41 * f(a + (6 * h))))
        A += 7

    return (S,A)

def make_plot(a, b, f, k, n, Ivalue, c, s, t):
    """
        Erzeugt Plots nach bestimmten Voraussetzungen

        Parameter:
        @param a:       Startpunkt des Intervalls (a,b)
        @type a:        int
        @param b:       Endpunkt des Intervalls (a,b)
        @type b:        int
        @param f:       die zu integrierende Funktion
        @param k:       k∈{1,2,6}, gibt an, welche k-te 
                        Newton–Cˆotes-Formel angewendet wird
                        Es gilt:
                        k=1: Trapezregel    1/2  1/2
                        k=2: Simpson-Regel  1/6  4/6  1/6
                        k=6: Weddle-Regel   41/840  216/840  27/840 
                                            272/840  27/840  216/840  41/840
        @type k:        int
        @param n:       Anzahl der Intervalle, auf denen 
                        die k-te Newton–Cˆotes-Formel 
                        (mit k∈{1,2,6}) angewendet wird
        @type n:        int
        @param Ivalue:  Wert für ein bestimmtes Integral
        @type Ivalue:   float
        @param c:       Farbe für den Grafen
        @type c:        string
        @param s:       Position für den subplot
        @type s:        int
        @param t:       Text für den Plottitel
        @type t:        string
    """

    x = list(np.arange(1,n+1,1))
    y = []
    y2 = []

    for i in x:
        # snc ist eine Liste von Tupeln (S, A) mit S = Quadraturwert und A = f-Auswertungen
        snc = summed_newton_cotes(a, b, f, k, i)
        y.append(abs(Ivalue - snc[0]))
        y2.append(snc[1])

    if s == 1:
        plt.subplot(2,1,s)
        plt.loglog(x, y, c, label='k=%s' % k) 
        plt.xlabel('n=%i (doppelte logarithmische Skala)' % n)
        plt.ylabel('Quadraturfehler')
        
    if s == 2:
        plt.subplot(2,1,s)
        # geplottet werden die f-Auswertungen über den Quadraturfehler
        # y2 = f-Auswertungen an der Stelle n (A(n))
        # y = Quadraturfehler an der Stelle n (S(n))
        # zur bessern Darstellung wird y2 auf der x-Achse abgebildet und
        # y auf der y-Achse
        plt.loglog(y2, y, c+'.', label='k=%s' % k)
        plt.xlabel('f-Auswertungen')
        plt.ylabel('Quadraturfehler')

    plt.tight_layout()
    plt.title(t)
    plt.legend()
    


# Variablen für Plot-window sin(x)
f1 = lambda x: math.sin(x)
I_f1_0_Pi = 2
a1 = 0
b1 = math.pi

# Variablen für Plot-window √x+ sin(21πx)
f2 = lambda x: math.sqrt(x) + math.sin(21 * math.pi * x)
I_f2_0_1 = (2 / 21) * (7 + (1 / math.pi))
a2 = 0
b2 = 1

n = 1000

plt.figure('Aufgabe 3b) Plots für sin(x)')
# b) Fehlerkurve mit n = 1, ..., 1000 für f(x) = sin(x)
# doppelt logarithmisch dargestellt
text = 'Fehlerkurve im Intervall [0, π].'
make_plot(a1, b1, f1, 1, n, I_f1_0_Pi, 'b', 1, text)
make_plot(a1, b1, f1, 2, n, I_f1_0_Pi, 'r', 1, text)
make_plot(a1, b1, f1, 6, n, I_f1_0_Pi, 'g', 1, text)

text = 'Quadrturfehler in Abhängigkeit der f-Auswertugen an der Stelle n'
make_plot(a1, b1, f1, 1, n, I_f1_0_Pi, 'b', 2, text)
make_plot(a1, b1, f1, 2, n, I_f1_0_Pi, 'r', 2, text)
make_plot(a1, b1, f1, 6, n, I_f1_0_Pi, 'g', 2, text)


plt.figure('Aufgabe 3c) Plots für √x+ sin(21πx)')
# c) Fehlerkurve mit n = 1, ..., 1000 für 
# doppelt logarithmisch dargestellt

text = 'Fehlerkurve im Intervall [0, 1].'
make_plot(a2, b2, f2, 1, n, I_f2_0_1, 'b', 1, text)
make_plot(a2, b2, f2, 2, n, I_f2_0_1, 'r', 1, text)
make_plot(a2, b2, f2, 6, n, I_f2_0_1, 'g', 1, text)

text = 'Quadrturfehler in Abhängigkeit der f-Auswertugen an der Stelle n'
make_plot(a2, b2, f2, 1, n, I_f2_0_1, 'b', 2, text)
make_plot(a2, b2, f2, 2, n, I_f2_0_1, 'r', 2, text)
make_plot(a2, b2, f2, 6, n, I_f2_0_1, 'g', 2, text)

plt.show()