from matplotlib import pyplot as plt

import math
import numpy as np

T = 10
x0 = 1
xQ0 = 1001

xvalues = np.linspace(0, T, 10*T)
yvalues01 = list(map(lambda t: math.exp(2 * t) * (1 + t**2), xvalues))
yvalues02 = list(map(lambda t: math.exp(2 * t) * (1001 + t**2), xvalues))

error = list(map(lambda x, y: abs(x - y), yvalues01, yvalues02))

plt.figure('Aufgabe 3b): e^(2t) · (x_0 + t^2)')
plt.subplot(2,1,1)
plt.loglog(xvalues, yvalues01, 'b', label='x_0 = 1')
plt.loglog(xvalues, yvalues02, 'r', label='x_0 = 1001')

plt.legend()
plt.xlabel('Anzahl x-Werte = %i (logarithmische Skala)' % (10*T))
plt.ylabel('f(x) \n(logarithmische Skala)')

plt.subplot(2,1,2)
plt.loglog(xvalues, error, label='|x(t) - x~(t)|')
plt.legend()
plt.xlabel('Anzahl x-Werte = %i (logarithmische Skala)' % (10*T))
plt.ylabel('errors \n(logarithmische Skala)')

plt.tight_layout()


plt.show()