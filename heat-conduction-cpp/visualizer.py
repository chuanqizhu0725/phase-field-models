import pandas as pd
import matplotlib.pyplot as plt

step_arr = [0, 200, 400]
for step in step_arr:
    df = pd.read_csv(f"data/1d{step}.csv", header=None)
    arr = df[0].values
    plt.plot(arr)
    plt.savefig(f"figures/1d{step}")
    plt.close()
