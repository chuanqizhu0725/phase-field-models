import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 64
ny = 32
df = pd.read_csv(f"data/con/2d40000.csv", header=None)
arr = df.values
mat = arr.reshape(nx, ny)
dft = pd.read_csv(f"data/temp/2d40000.csv", header=None)
arrt = dft[0].values
# plt.figure(figsize=(5, 5))
# plt.plot(mat[:, 8])
plt.subplot(1, 2, 1)
plt.imshow(mat)
plt.subplot(1, 2, 2)
plt.plot(arrt)
plt.show()
