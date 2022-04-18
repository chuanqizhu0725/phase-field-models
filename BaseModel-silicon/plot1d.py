import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 128
ny = 1
ns = 10
step_arr = np.arange(0, ns*10, ns)
for step in step_arr:
    df = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
    arr = df.values
    plt.plot(arr)
    plt.savefig(f"figures/phi/1d{step}.png")
    plt.close()

    # dfc = pd.read_csv(f"data/con/1d{step}.csv", header=None)
    # arrc = dfc.values
    # plt.plot(arrc)
    # plt.savefig(f"figures/con/1d{step}.png")
    # plt.close()

    dft = pd.read_csv(f"data/temp/1d{step}.csv", header=None)
    arrt = dft.values
    plt.plot(arrt)
    plt.savefig(f"figures/temp/1d{step}.png")
    plt.close()
