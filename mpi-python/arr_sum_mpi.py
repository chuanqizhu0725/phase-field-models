from mpi4py import MPI
import numpy as np
import time
t0 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

val = 0
arr = np.ones(12500000, "i")
for num in arr:
    val += num
val = np.array(val, "i")

value_sum = np.array(0, "i")
comm.Reduce(val, value_sum, op=MPI.SUM, root=0)

if rank == 0:
    print(value_sum)
    print(time.time()-t0)
