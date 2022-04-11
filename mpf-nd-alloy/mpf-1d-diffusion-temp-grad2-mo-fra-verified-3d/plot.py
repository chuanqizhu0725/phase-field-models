import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 2000
nx = 32
ny = 64
step_arr = np.arange(0, ns*101, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(nx, ny)
    plt.imshow(matc, origin='lower')
    plt.savefig(f"figures/con/2d{step}")
    plt.close()
