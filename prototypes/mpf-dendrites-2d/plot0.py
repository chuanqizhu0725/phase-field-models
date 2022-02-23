import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/con/2d20000.csv", header=None)
arr = df.values
mat = arr.reshape(200, 200)
plt.figure(figsize=(5, 5))
# plt.plot(mat[225, :])
plt.matshow(mat)
plt.show()
