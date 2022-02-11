###############################
#        Heat Conduction      #
# Finite Difference Method 2D #
###############################

# Boundary Condition: fixed value

import numpy as np
import matplotlib.pyplot as plt

Nx = 128
Ny = 128
dx = 1.0
dy = 1.0
dtime = 0.2
nstep = 600
nprint = 200

rad = 20

temp = np.zeros((Nx, Ny))

for i in range(Nx):
    for j in range(Ny):
        if (i-Nx/2)**2 + (j-Ny/2)**2 < rad**2:
            temp[i][j] = 1.0

for istep in range(nstep):
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            temp[i][j] = temp[i][j] + dtime * \
                (temp[i+1][j] + temp[i-1][j] + temp[i][j+1] +
                 temp[i][j-1] - 4.0*temp[i][j])/(dx*dx + dy*dy)
    if istep % nprint == 0:
        plt.matshow(temp)
        plt.show()
