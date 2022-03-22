import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ts = 40
step_arr = np.arange(0, ts*11, ts)
for step in step_arr:
    df = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
    arr = df[0].values
    plt.plot(arr)
    plt.savefig(f"figures/phi/1d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/con/1d{step}.csv", header=None)
    arrc = dfc[0].values
    # plt.plot(arr)
    plt.plot(arrc)
    plt.savefig(f"figures/con/1d{step}")
    plt.close()
