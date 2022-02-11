###############################
#        Heat Conduction      #
# Finite Difference Method 1D #
###############################

# Boundary Condition: fixed value

import numpy as np
import matplotlib.pyplot as plt
import time

t0 = time.time()

Nx = 128
dx = 1.0
dtime = 0.2
nstep = 5000
nprint = 200

temp = np.zeros(Nx)

for idx in range(Nx):
    if idx > 44 and idx < 84:
        temp[idx] = 1.0

for istep in range(nstep):
    for i in range(1, Nx-1):
        temp[i] = temp[i] + dtime*(temp[i+1] + temp[i-1] - 2.0*temp[i])/(dx*dx)
    # if istep % nprint == 0:
        # plt.plot(temp)
        # plt.show()

f = open("1d-ser.txt", "a")
for val in temp:
    f.write(str(val) + "\n")
f.close()

print(time.time()-t0)
