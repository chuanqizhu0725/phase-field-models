import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

step_arr = arr = np.arange(0, 10000, 1000)
for step in step_arr:
    df = pd.read_csv(f"data/test{step}.csv", header=None)
    arr = df.values
    mat = arr.reshape(400, 400)
    plt.figure(figsize=(5, 5))
    plt.matshow(mat)
    plt.savefig(f"figures/1d{step}")
    plt.close()
