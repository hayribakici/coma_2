from matplotlib import pyplot as plt
from scipy import special
from enum import IntEnum, Enum

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
        Löst ein AWP mittels des expliziten Euler-Verfahrens aus der VL:
        x_{k + 1} = x_k + τ(λx_k + f(t_k)),  mit k= 0,...,n - 1

        Parameter:
        @param A:   λ = A
        @type A:    int
        @param f:   Funktion fuer das Eulerverfahren: x'(t) = λx(t) + f(t)
        @param x0:  Startwert
        @type x0:   int ? float
        @param T:   Schrittweite τ := T / n
        @type T:    float
        @param n:   Anzahl der Euler-Auswertungen (Notwendig für T)
        @type n:    int

        @rtype:     [tupel] 
        @return:    Gibt eine Liste von Vektoren aller Iterierten als x ∈ ℝ^{n + 1} 
                    (inklusive Startwert) zurück
    """

    # Liste mit Startwert
    result = [(0,x0)]

    tau = _get_Schrittweite(T, n)

    for i in range(n):
        tk = result[i][0]
        xk = result[i][1]
        result.append(((tk + tau), (xk + tau * (A * xk + f(tk)))))

    return result

# TODO
# Zur Vorgehensweise: Die gegebene Funktion aus dem Skript (3.38) Seite 29 
# muss nach x_{k+1} aufgelöst werden. Für n=60 sieht der Näherungsgraf so aus, 
# wie auf dem Bild im Skript, aber für n=120 ist es sehr abenteuerlich und mit 
# loglog oder semilogy gar nicht richtig darzustellen. Ich nehme an, der Graf sieht ähnlich wie für n=60 aus...
# Aufgabenteil b)
def implicitEuler(A, f, x0, T, n):
    """
        Loest ein AWP mittels des impliziten Euler-Verfahrens aus der VL:
        x_{k + 1} = x_k + τ (λx_{k + 1} + f(t_{k + 1})),  mit k= 0,...,n - 1

        Um x_{k+1} in Abhängigkeit von x_{k+1} berechnen zu können, muss
        die Gleichung nach x_{k+1} aufgelöst werden

        => x_{k+1} = (x_k + tau ⋅ f(t_{k + 1})) / (1 - τ ⋅ λ)

        Parameter:
        @param A:   λ = A
        @type A:    int
        @param f:   Funktion über die Iteriert wird: x'(t) = λx(t) + f(t)
        @param x0:  Startwert
        @type x0:   int ? float
        @param T:   Schrittweite τ := T / n
        @type T:    float
        @param n:   Anzahl der Euler-Auswertungen (Notwendig für T)
        @type n:    int


        @rtype:     [tupel] 
        @return:    Gibt eine Liste von Vektoren aller Iterierten als x ∈ ℝ^{n + 1} 
                    (inklusive Startwert) zurueck in der Form (t,x) 
                    mit t = Wert auf x-Achse und x = Wert auf y-Achse
    """

    # Liste mit Startwert
    result = [(0,x0)]

    tau = _get_Schrittweite(T, n)

    for i in range(n):
        tk = result[i][0] + tau # == t_{k + 1}
        xk = result[i][1]
        result.append((tk, (xk + tau * f(tk)) / (1 - tau * A)))

    return result

# Aufgabenteil c)
def make_plot(A, f, x0, T, n, color, method_type, plot_position):
    """
        Erstellt Plots nach bestimmten Bedingungen

        Parameter:
        @param A:       λ = A
        @type A:        int
        @param f:       Funktion über die Iteriert wird: x'(t) = λx(t) + f(t)
        @param x0:      Startwert
        @type x0:       int ? float
        @param T:       Schrittweite τ := T / n
        @type T:        float
        @param n:       Anzahl der Euler-Auswertungen (Notwendig für T)
        @type n:        int
        @param c:       Farbe für den entsprechenden Grafen
        @type c:        string
        @param method:  EXPLICIT fuer explicit, IMPLICIT fuer implicit
        @type method:   MethodType
        @param sp:      Position des Plots
        @type sp:       PlotPosition
    """

    if method_type == MethodType.EXPLICIT:
        plotlist = explicitEuler(A, f, x0, T, n)

        plt.loglog(*zip(*plotlist), color, label='n = %i' % n)

        if n == 120:
            test = _get_test_fkt(testFkt, plotlist)
            plt.loglog(*zip(*test), '-r', label='Vergleichskurve n = %i' % n)

        plt.xlabel('Anzahl der x-Werte in Abhaengigkeit von n und der Schrittweite T/n (T = %i).' % T)
        plt.ylabel('Näherungswerte')
        plt.title('Explizites Eulerverfahren')

    elif method_type == MethodType.IMPLICIT:
        plotlist = implicitEuler(A, f, x0, T, n)

        if plot_position == PlotPosition.NONE:
            plotlist = implicitEuler(A, f, x0, T, n)
            plt.yscale('symlog')
            plt.semilogx(*zip(*plotlist), color, label='n = %i' % n)

            if n == 120:
                test = _get_test_fkt(testFkt, plotlist)
                plt.yscale('symlog')
                plt.semilogx(*zip(*test), '-r', label='Vergleichskurve n = %i' % n)

            plt.xlabel('Anzahl der x-Werte in Abhaengigkeit von n und der Schrittweite T/n (T = %i).' % T)
            plt.ylabel('Näherungswerte')
            plt.title('Implizites Eulerverfahren mit symlog der y-Achse und semilogx der x-Achse')

        if plot_position == PlotPosition.TOP:
            plt.subplot(3, 1, plot_position)
            plt.yscale('symlog')
            plt.semilogx(*zip(*plotlist), color, label='n = %i' % n)
            plt.title('Implizites Eulerverfahren mit n = %i' % n)
            plt.xlabel('Anzahl der x-Werte in Abhängigkeit von n und der Schrittweite T/n (T = %i).' % T)
            plt.ylabel('Näherungswerte')

        if plot_position == PlotPosition.CENTER:
            plt.subplot(3, 1, plot_position)
            plt.yscale('symlog')
            plt.semilogx(*zip(*plotlist), color, label='n = %i' % n)
            plt.title('Implizites Eulerverfahren mit n = %i' % n)
            plt.xlabel('Anzahl der x-Werte in Abhängigkeit von n und der Schrittweite T/n (T = %i).' % T)
            plt.ylabel('Näherungswerte')

            if n == 120:
                test = _get_test_fkt(testFkt, plotlist)
                plt.subplot(3, 1, plot_position)
                plt.yscale('symlog')
                plt.semilogx(*zip(*test), '-r', label='Vergleichskurve n = %i' % n)
                plt.title('Vergleichskurve mit n = %i' % n)

            plt.tight_layout()
    plt.legend()

class MethodType(Enum):
    IMPLICIT = 0,
    EXPLICIT = 1
        
# Using IntEnum to represent a restricted (from 0, ..., 3) int-type
class PlotPosition(IntEnum):
    """
    Position of the plot.
    """
    NONE = 0,
    TOP = 1,
    CENTER = 2,
    BOTTOM = 3


A = 16
f = lambda t: t
T = 8 
x0 = 1  # => Punkt (0, 1) als Startwert
testFkt = lambda t: (257 / 256) * math.exp(16 * t) - (t / 16) - (1 / 256)

# Plot fuer Explizites Eulerverfahren
plt.figure('Aufgabe 3b) Plots fuer das explizite Eulerverfahren mit f(t) = t')

# Explizites Eulerverfahren mit n = 60
n = 60
make_plot(A, f, x0, T, n, '-b', MethodType.EXPLICIT, PlotPosition.NONE)

# Explizites Eulerverfahren mit n = 120
n = 120
make_plot(A, f, x0, T, n, '-g', MethodType.EXPLICIT, PlotPosition.NONE)

# Plot fuer Implizites Eulerverfahren
plt.figure('Aufgabe 3b) Plots fuer das implizite Eulerverfahren mit f(t) = t')

# Implizites Eulerverfahren mit n = 60
n = 60
make_plot(A, f, x0, T, n, '-b', MethodType.IMPLICIT, PlotPosition.TOP)

# Implizites Eulerverfahren mit n = 120
n = 120
make_plot(A, f, x0, T, n, '-g', MethodType.IMPLICIT, PlotPosition.CENTER)

# Test eines gemeinsamen Plots fuer implizites Eulerverfahren
plt.figure('Aufgabe 3b) Plots fuer das implizite Eulerverfahren mit f(t) = t (Testplot)')

# Implizites Eulerverfahren mit n = 60
n = 60
make_plot(A, f, x0, T, n, '-b', MethodType.IMPLICIT, PlotPosition.NONE)

# Implizites Eulerverfahren mit n = 120
n = 120
make_plot(A, f, x0, T, n, '-g', MethodType.IMPLICIT, PlotPosition.NONE)

plt.show()