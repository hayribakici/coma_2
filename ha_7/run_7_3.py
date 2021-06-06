from matplotlib import pyplot as plt
from scipy import special

import math
import numpy as np

def _get_Schrittweite(T, n):
    return (T / n)

def _get_test_fkt(f, list):

    result = []

    for i in list:
        x = i[0]
        result.append((x, f(x)))

    return result

# Aufgabenteil a)
def explicitEuler(A, f, x0, T, n):
    """
        Loest ein AWP mittels des expliziten Euler-Verfahrens aus der VL:
        x_{k + 1} = x_k + tau(lambdax_k + f(t_k)),  mit k= 0,...,n - 1

        Parameter:
        @param A:   lambda = A
        @type A:    int
        @param f:   Funktion fuer das Eulerverfahren: x'(t) = lambdax(t) + f(t)
        @param x0:  Startwert
        @type x0:   int ? float
        @param T:   Schrittweite tau := T / n
        @type T:    float
        @param n:   Anzahl der Euler-Auswertungen (Notwendig fuer T)
        @type n:    int

        @rtype:     [tupel] 
        @return:    Gibt eine Liste von Vektoren aller Iterierten als x elem R^{n + 1} 
                    (inklusive Startwert) zurueck
    """

    # Liste mit Startwert
    result = [(0,x0)]

    h = _get_Schrittweite(T, n)

    for i in range(0, n):
        tk = result[i][0]
        xk = result[i][1]
        result.append(((tk + h), (xk + h * (A * xk + f(tk)))))

    return result

# TODO
# Aufgabenteil b)
def implicitEuler(A, f, x0, T, n):
    """
        Loest ein AWP mittels des impliziten Euler-Verfahrens aus der VL:
        x_{k + 1} = x_k + tau(lambdax_{k + 1} + f(t_{k + 1})),  mit k= 0,...,n - 1

        Um x_{k+1} in Abhaengigkeit von x_{k+1} berechnen zu koennen, muss
        die Gleichung nach x_{k+1} aufgeloest werden

        => x_{k+1} = (x_k + tau * f(t_{k + 1})) / (1 - tau * lambda)

        Parameter:
        @param A:   lambda = A
        @type A:    int
        @param f:   Funktion ueber die Iteriert wird: x'(t) = lambdax(t) + f(t)
        @param x0:  Startwert
        @type x0:   int ? float
        @param T:   Schrittweite tau := T / n
        @type T:    float
        @param n:   Anzahl der Euler-Auswertungen (Notwendig fuer T)
        @type n:    int


        @rtype:     [tupel] 
        @return:    Gibt eine Liste von Vektoren aller Iterierten als x elem R^{n + 1} 
                    (inklusive Startwert) zurueck in der Form (t,x) 
                    mit t = Wert auf x-Achse und x = Wert auf y-Achse
    """

    # Liste mit Startwert
    result = [(0,x0)]

    h = _get_Schrittweite(T, n)

    for i in range(0, n):
        tk = result[i][0] + h
        xk = result[i][1]
        result.append(((tk, ((xk + h * f(tk)) / (1 - h * A)))))

    return result

# Aufgabenteil c)
def make_plot(A, f, x0, T, n, c, method, sp):
    """
        Erstellt Plots nach bestimmten Bedingungen

        Parameter:
        @param A:       lambda = A
        @type A:        int
        @param f:       Funktion ueber die Iteriert wird: x'(t) = lambdax(t) + f(t)
        @param x0:      Startwert
        @type x0:       int ? float
        @param T:       Schrittweite tau := T / n
        @type T:        float
        @param n:       Anzahl der Euler-Auswertungen (Notwendig fuer T)
        @type n:        int
        @param c:       Farbe fuer den entsprechenden Grafen
        @type c:        string
        @param method:  e fuer explicit, i fuer implicit
        @type method:   char bzw. string
        @param sp:      Position fuer den subplot: 0 = kein subplot, 
                        1 = oben, 2 = mitte, 3 = unten
        @type sp:       int
    """

    if method == 'e':
        plotlist = explicitEuler(A, f, x0, T, n)

        plt.loglog(zip(plotlist), c, label='n = %i' % n)

        if n == 120:
            test = _get_test_fkt(testFkt, plotlist)
            plt.loglog(zip(test), '-r', label='Vergleichskurve n = %i' % n)

        plt.xlabel('Anzahl der x-Werte in Abhaengigkeit von n und der Schrittweite T/n (T = %i).' % T)
        plt.ylabel('Naeherungswerte')
        plt.title('Explizites Eulerverfahren')
        plt.legend()

    if sp == 0 and method == 'i':
        plotlist = implicitEuler(A, f, x0, T, n)
        plt.yscale('symlog')
        plt.semilogx(zip(plotlist), c, label='n = %i' % n)

        if n == 120:
            test = _get_test_fkt(testFkt, plotlist)
            plt.yscale('symlog')
            plt.semilogx(*zip(*test), '-r', label='Vergleichskurve n = %i' % n)

        plt.xlabel('Anzahl der x-Werte in Abhaengigkeit von n und der Schrittweite T/n (T = %i).' % T)
        plt.ylabel('Naeherungswerte')
        plt.title('Implizites Eulerverfahren mit symlog der y-Achse und semilogx der x-Achse')
        plt.legend()
    
    if sp > 0 and method == 'i':
        plotlist = implicitEuler(A, f, x0, T, n)
        if sp == 1:
            plt.subplot(3,1,sp)
            plt.yscale('symlog')
            plt.semilogx(*zip(*plotlist), c, label='n = %i' % n)
            plt.title('Implizites Eulerverfahren mit n = %i' % n)
            plt.legend()
            plt.xlabel('Anzahl der x-Werte in Abhaengigkeit von n und der Schrittweite T/n (T = %i).' % T)
            plt.ylabel('Naeherungswerte')

        if sp == 2:
            plt.subplot(3,1,sp)
            plt.yscale('symlog')
            plt.semilogx(*zip(*plotlist), c, label='n = %i' % n)
            plt.title('Implizites Eulerverfahren mit n = %i' % n)
            plt.legend()
            plt.xlabel('Anzahl der x-Werte in Abhaengigkeit von n und der Schrittweite T/n (T = %i).' % T)
            plt.ylabel('Naeherungswerte')

        if n == 120:
            test = _get_test_fkt(testFkt, plotlist)
            plt.subplot(3,1,3)
            plt.yscale('symlog')
            plt.semilogx(*zip(*test), '-r', label='Vergleichskurve n = %i' % n)
            plt.title('Vergleichskurve mit n = %i' % n)
            plt.legend()

        plt.tight_layout()

        
        

A = 16
f = lambda t: t
T = 8 
x0 = 1  # => Punkt (0, 1) als Startwert
testFkt = lambda t: (257 / 256) * math.exp(16 * t) - (t / 16) - (1 / 256)

# Plot fuer Explizites Eulerverfahren
plt.figure('Aufgabe 3b) Plots fuer das explizite Eulerverfahren mit f(t) = t')

# Explizites Eulerverfahren mit n = 60
n = 60
make_plot(A, f, x0, T, n, '-b', 'e', 0)

# Explizites Eulerverfahren mit n = 120
n = 120
make_plot(A, f, x0, T, n, '-g', 'e', 0)

# Plot fuer Implizites Eulerverfahren
plt.figure('Aufgabe 3b) Plots fuer das implizite Eulerverfahren mit f(t) = t')

# Implizites Eulerverfahren mit n = 60
n = 60
make_plot(A, f, x0, T, n, '-b', 'i', 1)

# Implizites Eulerverfahren mit n = 120
n = 120
make_plot(A, f, x0, T, n, '-g', 'i', 2)

# Test eines gemeinsamen Plots fuer implizites Eulerverfahren
plt.figure('Aufgabe 3b) Plots fuer das implizite Eulerverfahren mit f(t) = t (Testplot)')

# Implizites Eulerverfahren mit n = 60
n = 60
make_plot(A, f, x0, T, n, '-b', 'i', 0)

# Implizites Eulerverfahren mit n = 120
n = 120
make_plot(A, f, x0, T, n, '-g', 'i', 0)

plt.show()