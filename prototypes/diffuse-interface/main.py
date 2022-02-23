import numpy as np
import matplotlib.pyplot as plt
import time

Nx = 128
dx = 1.0
dtime = 0.2
nstep = 5000
nprint = 200

phi = np.zeros(Nx)

for idx in range(Nx):
    if idx > 44 and idx < 84:
        phi[idx] = 1.0

for istep in range(nstep):
    for i in range(1, Nx-1):
        phi[i] = phi[i] + dtime*(phi[i+1] + phi[i-1] - 2.0*phi[i])/(dx*dx)
    if istep % nprint == 0:
        plt.plot(phi)
        plt.show()
