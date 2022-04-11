import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/con/2d200.csv", header=None)
arr = df.values
mat = arr.reshape(64, 64)
plt.figure(figsize=(5, 5))
# plt.plot(mat[:, 50])
plt.matshow(mat)
plt.show()
