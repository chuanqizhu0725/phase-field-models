import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/con/2d0.csv", header=None)
arr = df.values
mat = arr.reshape(128, 128)
# plt.figure(figsize=(5, 5))
plt.plot(mat[:, 110])
# plt.matshow(mat)
plt.show()
