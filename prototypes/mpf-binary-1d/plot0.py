import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/con/2d500.csv", header=None)
arr = df.values
mat = arr.reshape(100, 100)
plt.figure(figsize=(10, 2))
plt.plot(mat[25, :])
# plt.matshow(mat)
plt.show()
