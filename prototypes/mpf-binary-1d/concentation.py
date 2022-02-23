import numpy as np
import matplotlib.pyplot as plt

ch = np.zeros(100)

cah = np.zeros(100)
cbh = np.zeros(100)

sah = np.zeros(100)
sbh = np.zeros(100)

for i in range(100):
    if i < 50:
        cah[i] = 0.1
    elif i > 49:
        cbh[i] = 0.2

# plt.plot(sah)
plt.plot(cah)
plt.plot(cbh)
# plt.plot(ch)
plt.show()


# x = 1.0-(i-44)*0.1
# sah[i] = x**3*(10.0-15.0*x+6.0*x**2)
