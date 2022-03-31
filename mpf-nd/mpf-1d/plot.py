import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 40
step_arr = np.arange(0, ns*11, ns)
for step in step_arr:
    df = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
    arr = df[0].values
    plt.plot(arr)
    plt.savefig(f"figures/phi/1d{step}")
    plt.close()
