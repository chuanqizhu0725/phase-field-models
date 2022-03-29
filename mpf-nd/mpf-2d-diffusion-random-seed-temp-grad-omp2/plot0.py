import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/con/2d4800.csv", header=None)
arr = df.values
mat = arr.reshape(128, 128)
plt.figure(figsize=(5, 5))
# plt.plot(mat[36, :])
plt.matshow(mat)
plt.show()
