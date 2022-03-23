import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 64
df = pd.read_csv(f"data/con/2d400.csv", header=None)
arr = df.values
mat = arr.reshape(nx, nx)
# plt.plot(mat[:, 32])
plt.figure(figsize=(5, 5))
plt.matshow(mat)
plt.show()
