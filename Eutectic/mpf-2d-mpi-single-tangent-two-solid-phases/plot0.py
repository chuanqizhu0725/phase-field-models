import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("phi0.dat", header=None)
arr = df.values
mat = arr.reshape(100, 100)
plt.figure(figsize=(5, 5))
# plt.plot(mat[:, 200])
plt.matshow(mat)
plt.show()
