import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

step_arr = np.arange(0, 110000, 10000)
for step in step_arr:
    df = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
    arr = df[0].values
    plt.plot(arr)
    plt.savefig(f"figures/phi/1d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/con/1d{step}.csv", header=None)
    arrc = dfc[0].values
    dfcs = pd.read_csv(f"data/cons/1d{step}.csv", header=None)
    arrcs = dfcs[0].values
    dfcl = pd.read_csv(f"data/conl/1d{step}.csv", header=None)
    arrcl = dfcl[0].values
    plt.plot(arr)
    plt.plot(arrc)
    # plt.plot(arrcs)
    # plt.plot(arrcl)
    plt.savefig(f"figures/con/1d{step}")
    plt.close()
