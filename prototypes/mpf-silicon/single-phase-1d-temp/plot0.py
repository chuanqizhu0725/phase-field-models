import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/temp/1d900000.csv", header=None)
arr = df.values
# mat = arr.reshape(400, 400)
# plt.figure(figsize=(5, 5))
# plt.plot(mat[225, :])
# plt.matshow(mat)
plt.plot(arr)
plt.show()
