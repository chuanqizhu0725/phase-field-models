import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 50000
nx = 80
ny = 80
step_arr = np.arange(0, ns*11, ns)
for step in step_arr:
    df = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
    arr = df[0].values
    mat = arr.reshape(nx, ny)
    plt.imshow(mat)
    plt.savefig(f"figures/phi/2d{step}")
    plt.close()

    dft = pd.read_csv(f"data/temp/2d{step}.csv", header=None)
    arrt = dft[0].values
    matt = arrt.reshape(nx, ny)
    plt.imshow(matt)
    plt.savefig(f"figures/temp/2d{step}")
    plt.close()
