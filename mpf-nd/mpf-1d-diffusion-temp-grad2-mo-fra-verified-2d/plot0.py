import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = ny = 128
df = pd.read_csv(f"data/con/2d1000.csv", header=None)
arr = df.values
mat = arr.reshape(nx, ny)
# plt.figure(figsize=(5, 5))
plt.plot(mat[:, 64])
# plt.matshow(mat)
plt.show()
