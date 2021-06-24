import numpy as np

# 2b
def taylor2(mu, x0, T, tau):

    n = 10
    return __taylor2(mu, x0, T, tau, n - 1)

def __taylor2(mu, x0, T, tau, k):
    
    if k <= 0:
        return x0
    
    return __taylor2(mu, x0, T, tau, k - 1) - (tau * mu + 1) * __taylor2(mu, x0, T, tau, k - 2)
    


mus = [-1, -2, -3]
T = 1
x0 = 1
tau = 999999999999

for mu in mus:
    print("τ = ", tau, "μ =", mu, ":", taylor2(mu, x0, T, tau))