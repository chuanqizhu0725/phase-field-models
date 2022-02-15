###############################
#     Equillibrium shape     #
# Finite Difference Method 2D #
###############################

# Boundary Condition: Periodic

import numpy as np
import matplotlib.pyplot as plt
import math

# parameters of calculation
Nx = 100
Ny = 100
NxNy = Nx*Ny
dx = 0.03
dy = 0.03
dtime = 1.0e-4
nstep = 4000
nprint = 1000

# parameters of materials
tau = 0.0003
epsilonb = 0.01
delta = 0.02
aniso = 4.0
theta0 = math.pi/4.0

# initialize the seed
phi = np.zeros((Nx, Ny))
lap_phi = np.zeros((Nx, Ny))
phidx = np.zeros((Nx, Ny))
phidy = np.zeros((Nx, Ny))
epsilon = np.zeros((Nx, Ny))
epsilon_deriv = np.zeros((Nx, Ny))
rad = 15

for i in range(Nx):
    for j in range(Ny):
        if (i-Nx/2)**2 + (j-Ny/2)**2 < rad**2:
            phi[i][j] = 1.0

for istep in range(nstep):
    if (istep == 0) or ((istep+1) % nprint == 0):
        plt.matshow(phi)
        plt.savefig(f"figure-delta002-m005/2d{istep+1}")
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
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

            # laplacian of phi
            lap_phi[i][j] = (phi[im][j] + phi[ip][j] + phi[i]
                             [jm] + phi[i][jp] - 4.0*phi[i][j])/(dx*dy)

            # gradients of phi
            phidx[i][j] = (phi[ip][j] - phi[im][j])/(2.0*dx)
            phidy[i][j] = (phi[i][jp] - phi[i][jm])/(2.0*dx)

            if phidx[i][j] == 0.0:
                theta = math.atan(phidy[i][j]/1.0e-6)
            elif phidx[i][j] == 0.0 and phidy[i][j] > 0.0:
                theta = math.pi/2.0
            elif phidx[i][j] == 0.0 and phidy[i][j] < 0.0:
                theta = -1.0*math.pi/2
            elif phidx[i][j] == 0.0 and phidy[i][j] == 0.0:
                theta = 0.0

            # epsilon and its derivative
            epsilon[i][j] = epsilonb * \
                (1.0 + delta*math.cos(aniso*(theta-theta0)))
            epsilon_deriv[i][j] = -1*epsilonb*aniso * \
                delta*math.sin(aniso*(theta-theta0))

    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
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

            term1 = (epsilon[i][jp]*epsilon_deriv[i][jp]*phidx[i][jp] -
                     epsilon[i][jm]*epsilon_deriv[i][jm]*phidx[i][jm])/(2.0*dy)

            term2 = (epsilon[ip][j]*epsilon_deriv[ip][j]*phidx[ip][j] -
                     epsilon[im][j]*epsilon_deriv[im][j]*phidx[im][j])/(2.0*dx)

            phi[i][j] = phi[i][j] + (dtime/tau)*(term1 + term2 + (
                epsilon[i][j]**2)*lap_phi[i][j] + phi[i][j]*(1.0-phi[i][j])*(phi[i][j]-0.45))
