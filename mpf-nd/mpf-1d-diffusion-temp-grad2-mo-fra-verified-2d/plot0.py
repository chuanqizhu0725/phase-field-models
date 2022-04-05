import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 128
ny = 64
df = pd.read_csv(f"data/con/2d50000.csv", header=None)
arr = df.values
mat = arr.reshape(nx, ny)
dft = pd.read_csv(f"data/temp/2d50000.csv", header=None)
arrt = dft[0].values
# plt.figure(figsize=(5, 5))
# plt.plot(mat[:, 8])
plt.subplot(1, 2, 1)
plt.imshow(mat, origin='lower')
plt.subplot(1, 2, 2)
plt.plot(arrt)
plt.show()
