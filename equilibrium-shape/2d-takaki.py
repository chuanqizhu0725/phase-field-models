###############################
#     Equillibrium shape     #
# Finite Difference Method 2D #
###############################

# Boundary Condition: Periodic

import numpy as np
import matplotlib.pyplot as plt
import math
import time

t0 = time.time()

# parameters of calculation
Nx = 64
Ny = 64
NxNy = Nx*Ny
dx = 0.3
dy = 0.3

dtime = 0.01
nstep = 1000
nprint = 200

# parameters of materials
delta = dx*5.0
gamma = 0.37e7
bb = 2.2
astre = 0.03

epsilonb = math.sqrt(3.0*delta*gamma/bb)
www = 6.0*gamma*bb/delta
mobi = 3.0

print(epsilonb)

phi = np.zeros((Nx, Ny), dtype=np.float64)
lap_phi = np.zeros((Nx, Ny), dtype=np.float64)
epsilon = np.zeros((Nx, Ny), dtype=np.float64)
phidx = np.zeros((Nx, Ny), dtype=np.float64)
phidy = np.zeros((Nx, Ny), dtype=np.float64)
termx = np.zeros((Nx, Ny), dtype=np.float64)
termy = np.zeros((Nx, Ny), dtype=np.float64)
rad = 15

for i in range(Nx):
    for j in range(Ny):
        if (i-Nx/2)**2 + (j-Ny/2)**2 < rad**2:
            phi[i][j] = 1.0

for istep in range(nstep):
    if istep == 0:
        print("*****************************\ncomputation starts!")
    if istep % 200 == 0 and istep > 0:
        print(
            f"run {istep-199}~{istep} steps of {nstep}, {time.time()-t0} has passed")
    if (istep == 0) or ((istep+1) % nprint == 0):
        plt.matshow(phi)
        plt.colorbar()
        plt.savefig(f"figure-takaki/2d{istep+1}")
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
            phidx = (phi[ip][j] - phi[im][j])/(2.0*dx)
            phidy = (phi[i][jp] - phi[i][jm])/(2.0*dy)

            if lap_phi[i][j] != 0.0 and phidx+phidy != 0.0:

                # gradient coefficient
                epsilon[i][j] = epsilonb*(1.0-3.0*astre)*(1.0+4.0*astre /
                                                          (1.0-3.0*astre)*(phidx**4+phidy**4)/(phidx+phidy)**4)

                # termx[i][j] = epsilon[i][j]*4.0*astre*(4*(phidx**3) *
                #                                        abs(phidx+phidy)-4.0*(phidx**4 + phidy**4))/abs(phidx+phidy)**5*(phidx+phidy)**2

                # termy[i][j] = epsilon[i][j]*4.0*astre*(4*(phidy**3) *
                #                                        abs(phidx+phidy)-4.0*(phidx**4 + phidy**4))/abs(phidx+phidy)**5*(phidx+phidy)**2

    plt.matshow(epsilon)
    plt.colorbar()
    plt.show()

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

            if lap_phi[i][j] != 0.0:

                term1 = (termx[ip][j]-termx[im][j])/(2.0*dx)
                term2 = (termy[i][jp]-termy[i][jm])/(2.0*dy)

                phi[i][j] = phi[i][j] + dtime*mobi*(epsilon[i][j]**2*lap_phi[i][j] + term1 + term2 +
                                                    4*www*phi[i][j]*(1.0-phi[i][j])*(phi[i][j]-0.45))
