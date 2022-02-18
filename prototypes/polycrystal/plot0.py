import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv(f"data/test5000.csv", header=None)
arr = df.values
mat = arr.reshape(400, 400)
plt.figure(figsize=(5, 5))
plt.matshow(mat)
plt.show()
