# https://rabernat.github.io/research_computing/parallel-programming-with-mpi-for-python.html
from mpi4py import MPI
import numpy as np

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# print("My rank is", rank)

# if rank == 0:
#     data = "hello, I am process 0!"
#     comm.send(data, dest=1)
# elif rank == 1:
#     data = comm.recv(source=0)
#     print("The process 1 received a msg", data)

# Create some np arrays on each process:
# For this demo, the arrays have only one
# entry that is assigned to be the rank of the processor
# value = np.array(rank, 'd')

# print(' Rank: ', rank, ' value = ', value)

# # initialize the np arrays that will store the results:
# value_sum = np.array(0.0, 'd')
# value_max = np.array(0.0, 'd')

# # perform the reductions:
# comm.Reduce(value, value_sum, op=MPI.SUM, root=0)
# comm.Reduce(value, value_max, op=MPI.MAX, root=0)

# if rank == 0:
#     print(' Rank 0: value_sum =    ', value_sum)
#     print(' Rank 0: value_max =    ', value_max)

# if rank == 0:
#     data = "hello"
# else:
#     data = None

# data = comm.bcast(data, root=0)
# print('Rank: ', rank, ', data: ', data)

# if rank == 0:
#     # create a data array on process 0
#     # in real code, this section might
#     # read in data parameters from a file
#     numData = 10
#     data = np.linspace(0.0, 3.14, numData)
# else:
#     numData = None

# # broadcast numData and allocate array on other ranks:
# numData = comm.bcast(numData, root=0)
# if rank != 0:
#     data = np.empty(numData, dtype='d')

# comm.Bcast(data, root=0)  # broadcast the array from rank 0 to all others

# print('Rank: ', rank, ', data received: ', data)

globalData = None

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i in range(2):
    numDataPerRank = 10
    sendbuf = np.linspace(rank*numDataPerRank+1, (rank+1)
                          * numDataPerRank, numDataPerRank)
    print('Rank: ', rank, ', sendbuf: ', sendbuf)

    recvbuf = None
    if rank == 0:
        recvbuf = np.empty(numDataPerRank*size, dtype='d')

    comm.Gather(sendbuf, recvbuf, root=0)

    if rank == 0:
        # globalData = recvbuf
        print('Rank: ', rank, ', recvbuf received: ', recvbuf)

print(globalData)
