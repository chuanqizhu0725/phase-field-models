import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv("undercooling.csv", header=None, sep=' ')
sp = df.loc[0].values
v1 = df.loc[1].values
v2 = df.loc[2].values
v3 = df.loc[3].values
fig, ax = plt.subplots(1, 3, figsize=(9, 3))
fig.suptitle(
    'Interface undercooling for varying lamellar spacing under different velocities', fontsize=16)
fig.supxlabel("(*) nondimensional value", x=0.8, y=0.015, fontsize=10)
ax[0].plot(sp, v1, marker='+')
ax[0].set_ylabel("undercooling", fontsize=14)
ax[0].set_xlabel("Lamellar spacing in grid number", fontsize=14)
ax[0].set_title("v=1.0$^{(*)}$", x=0.5, y=0.9)
ax[1].plot(sp, v2, marker='+')
ax[1].set_title("v=0.75", x=0.5, y=0.9)
ax[2].plot(sp, v3, marker='+')
ax[2].set_title("v=0.5", x=0.5, y=0.9)
plt.show()
