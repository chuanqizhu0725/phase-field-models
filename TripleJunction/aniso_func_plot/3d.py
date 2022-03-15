from matplotlib import pyplot as plt
import numpy as np


ep = 0.2

ax = plt.axes(projection='3d')

theta = np.linspace(0, 2*np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)

X = np.zeros(10000)
Y = np.zeros(10000)
Z = np.zeros(10000)
R = np.zeros(10000)

i = 0
for t in theta:
    for p in phi:
        x0 = np.sin(p)*np.cos(t)
        y0 = np.sin(p)*np.sin(t)
        z0 = np.cos(p)
        r = 1-3*ep+4*ep * \
            (np.sin(p)**4*(np.sin(t)**4 + np.cos(t)**4)+np.cos(p)**4)
        X[i] = x0*r
        Y[i] = y0*r
        Z[i] = z0*r
        R[i] = r
        i += 1

ax.scatter(X, Y, Z, c=R, cmap='viridis', linewidth=0.5)

plt.show()
