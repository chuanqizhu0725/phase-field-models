import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 128
ny = 56
# nxarr = np.arange(0, 128, 1)
df = pd.read_csv(f"data/con/2d50000.csv", header=None)
arr = df.values
mat = arr.reshape(nx, ny)
dft = pd.read_csv(f"data/temp/2d50000.csv", header=None)
arrt = dft[0].values
fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]})
ax[0].plot(arrt)
ax[0].set_ylabel("Temperature")
ax[0].set_xlabel("Height / Grid number")
ax[1].imshow(mat, origin='lower')
ax[1].set_ylabel("Height / Grid number")
ax[1].set_xlabel("Width / Grid number")
plt.show()
