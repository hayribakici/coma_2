from matplotlib import pyplot as plt

import math
import numpy as np

T = 10
# zum Zeichnen einer schoeneren
T_ = 100*T
# Startwerte: 
# x_0 = 1 und x~_0 = 1001

# Generieren der x-Werte
xvalues = np.linspace(0, T, T_)
# Berechnung der y-Werte fuer x_0 und x~_0 fuer f(x) = e^(2t) * (x_0 + t^2)
yvalues01 = list(map(lambda t: math.exp(2 * t) * (1 + t**2), xvalues))
yvalues02 = list(map(lambda t: math.exp(2 * t) * (1001 + t**2), xvalues))

# Berechnung der y-Werte fuer den Fehler |x(t) - x~(t)| 
error = list(map(lambda x, y: abs(x - y), yvalues01, yvalues02))

plt.figure('Aufgabe 3b): e^(2t) * (x_0 + t^2)')
plt.subplot(2,1,1)
# Grafen fuer x und x~
plt.loglog(xvalues, yvalues01, 'b', label='x_0 = 1')
plt.loglog(xvalues, yvalues02, 'r', label='x~_0 = 1001')

plt.legend()
plt.xlabel('Anzahl t-Werte = %i (logarithmische Skala)' % T_)
plt.ylabel('x(t) \n(logarithmische Skala)')

plt.subplot(2,1,2)
# Graf fuer |x(t) - x~(t)|
plt.loglog(xvalues, error, label='|x(t) - x~(t)|')

plt.legend()
plt.xlabel('Anzahl t-Werte = %i (logarithmische Skala)' % T_)
plt.ylabel('errors \n(logarithmische Skala)')

plt.tight_layout()

plt.show()