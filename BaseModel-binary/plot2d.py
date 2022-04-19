import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 128
ny = 128
ns = 100
step_arr = np.arange(0, ns*10, ns)
for step in step_arr:
    df = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
    arr = df[0].values
    mat = arr.reshape(nx, ny)
    plt.imshow(mat)
    plt.savefig(f"figures/phi/2d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(nx, ny)
    plt.imshow(matc)
    plt.savefig(f"figures/con/2d{step}")
    plt.close()
