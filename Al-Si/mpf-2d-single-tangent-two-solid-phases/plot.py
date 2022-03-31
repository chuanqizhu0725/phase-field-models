import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 64
step_arr = np.arange(0, 1100, 100)
for step in step_arr:
    df = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arr = df[0].values
    mat = arr.reshape(nx, nx)
    plt.figure(figsize=(5, 5))
    plt.matshow(mat)
    plt.savefig(f"figures/con/2d{step}")
    plt.close()
    # plt.show()
