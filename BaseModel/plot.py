import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 128
ny = 1
df = pd.read_csv("2d.csv", header=None)
arr = df.values
mat = arr.reshape(nx, ny)
# plt.figure(figsize=(5, 5))
# plt.matshow(mat)
plt.plot(arr)
plt.show()
