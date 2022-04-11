import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 2000
nx = 128
ny = 64
step_arr = np.arange(0, ns*101, ns)
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
    dft = pd.read_csv(f"data/temp/2d{step}.csv", header=None)
    arrt = dft[0].values
    plt.subplot(1, 2, 1)
    plt.imshow(matc, origin='lower')
    plt.subplot(1, 2, 2)
    plt.plot(arrt)
    plt.suptitle(f"{step} steps")
    plt.savefig(f"figures/con/2d{step}")
    plt.close()
