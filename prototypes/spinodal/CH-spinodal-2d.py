###############################
#    Spinodal Decomposition   #
# Finite Difference Method 2D #
###############################

# Boundary Condition: periodic

import numpy as np
import matplotlib.pyplot as plt
import random

Nx = 64
Ny = 64
dx = 1.0
dy = 1.0
dtime = 0.02
nstep = 1000
nprint = 250

AA = 1.0
c0 = 0.40
mob = 1.0
grad_coef = 0.5

con = np.zeros((Nx, Ny))
dummy = np.zeros((Nx, Ny))

for i in range(Nx):
    for j in range(Ny):
        con[i][j] = c0 + 0.02*(0.5-random.random())

plt.matshow(con)
plt.savefig("ch-spinodal-2d-img/0000.png")

for istep in range(nstep):
    for i in range(Nx):
        for j in range(Ny):
            jp = j+1
            jm = j-1
            ip = i+1
            im = i-1
            if jp == Ny:
                jp = 0
            if jm == -1:
                jm = Ny-1
            if ip == Nx:
                ip = 0
            if im == -1:
                im = Nx-1
            lap_con_ij = (con[im][j] + con[ip][j] + con[i]
                          [jm] + con[i][jp] - 4.0*con[i][j])/(dx*dy)
            f_con_ij = 2.0*AA*con[i][j]*(1-con[i][j])*(1-2.0*con[i][j])
            dummy[i][j] = f_con_ij - grad_coef*lap_con_ij

    for i in range(Nx):
        for j in range(Ny):
            jp = j+1
            jm = j-1
            ip = i+1
            im = i-1
            if jp == Ny:
                jp = 0
            if jm == -1:
                jm = Ny-1
            if ip == Nx:
                ip = 0
            if im == -1:
                im = Nx-1
            con[i][j] = con[i][j] + dtime*mob * \
                (dummy[im][j]+dummy[ip][j]+dummy[i][jm] +
                 dummy[i][jp]-4.0*dummy[i][j])/(dx*dy)
            if con[i][j] > 0.99999:
                con[i][j] = 0.99999
            if con[i][j] < 0.00001:
                con[i][j] = 0.00001
    if (istep+1) % nprint == 0 and istep > 0:
        plt.matshow(con)
        plt.savefig(f"ch-spinodal-2d-img/{istep+1}.png")
