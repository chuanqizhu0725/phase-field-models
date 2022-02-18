import numpy as np
import time

t0 = time.time()
arr = np.ones(100000000, "i")
sum = 0
for num in arr:
    sum += num
print(sum)
print(time.time()-t0)
