import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 32
ny = 32
# nxarr = np.arange(0, 128, 1)
df = pd.read_csv(f"data/con2d/2d2000.csv", header=None)
arr = df.values
mat = arr.reshape(nx, ny)
plt.imshow(mat)
plt.show()
