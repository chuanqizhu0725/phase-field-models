import numpy as np
import matplotlib.pyplot as plt
import math

epsilonb = 0.016
delta = 0.15
aniso = 4.0
theta0 = math.pi/4.0

theta = np.arange(0, 2*math.pi, 0.01)
r = np.zeros(len(theta))

for idx, val in enumerate(theta):
    r[idx] = epsilonb*(1.0 + delta*math.cos(aniso*(val-theta0)))

plt.subplots(subplot_kw={'projection': 'polar'})
plt.plot(theta, r)
plt.show()
