import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/temp/2d35000.csv", header=None)
arr = df.values
mat = arr.reshape(250, 250)
# plt.figure(figsize=(10, 2))
# plt.plot(mat[200, :])
plt.figure(figsize=(5, 5))
plt.matshow(mat)
plt.show()
