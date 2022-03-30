import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/con/1d12800.csv", header=None)
arr = df.values
df0 = pd.read_csv(f"data/con/1d12800.txt", header=None)
arr0 = df0.values
# plt.plot(arr)
plt.plot(arr0)
plt.show()
