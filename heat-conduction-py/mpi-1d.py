###############################
#        Heat Conduction      #
# Finite Difference Method 1D #
###############################

# Boundary Condition: fixed value
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import time

t0 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

Nx = 128
Nx2 = 64
dx = 1.0
dtime = 0.2
nstep = 5000
nprint = 200

temp = np.zeros(Nx)
temp00 = np.zeros(Nx2)
temp11 = np.zeros(Nx2)

for idx in range(Nx):
    if idx > 44 and idx < 84:
        temp[idx] = 1.0

for istep in range(nstep):
    if rank == 0:
        # the index of the rightmost grid is Nx/2-1
        if(istep == 0):
            temp0 = temp[0:Nx2]
        else:
            temp0 = temp00
        dataSend = temp0[Nx2-1]
        comm.send(dataSend, dest=1)
        dataRecv = comm.recv(source=1)
        for i in range(1, Nx2):
            west = temp0[i-1]
            center = temp0[i]
            east = 0.0
            if i == Nx2 - 1:
                east = dataRecv
            else:
                east = temp0[i+1]
            temp00[i] = center + dtime*(east + west - 2.0*center)/(dx*dx)

    if rank == 1:
        # the index of the leftmost grid is Nx/2
        if(istep == 0):
            temp1 = temp[Nx2:Nx]
        else:
            temp1 = temp11
        sendBuf = np.zeros(Nx2)
        dataSend = temp1[0]
        comm.send(dataSend, dest=0)
        dataRecv = comm.recv(source=0)
        for i in range(0, Nx2-1):
            east = temp1[i+1]
            center = temp1[i]
            west = 0.0
            if i == 0:
                west = dataRecv
            else:
                west = temp1[i-1]
            temp11[i] = center + dtime*(east + west - 2.0*center)/(dx*dx)

recvBuf = None
if rank == 0:
    sendBuf = temp00
    recvBuf = np.empty(Nx, dtype='d')

if rank == 1:
    sendBuf = temp11

comm.Gather(sendBuf, recvBuf, root=0)

if rank == 0:
    f = open("1d-mpi.txt", "a")
    for val in recvBuf:
        f.write(str(val) + "\n")
    f.close()

print(time.time()-t0)

# plt.plot(recvBuf)
# plt.show()
