#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define NX 20
#define NY 20
#define STEPS 100
#define BEGIN 1
#define UTAG 2
#define DTAG 3
#define NONE 0
#define DONE 4
#define MASTER 0

struct Parms
{
    float cx;
    float cy;
} parms = {0.1, 0.1};

int main(int argc, char *argv[])
{
    void inidat(), update(), prtdat();

    int taskid,
        numworkers,
        numtasks,
        rows, offset,
        dest, source,
        up, down,
        msgtype,
        rc, start, end,
        i, ix, iy, iz, it;

    MPI_Status status;

    // Allocate taskid to each core (cores = tasks = master + workers)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    numworkers = numtasks - 1;

    rows = NX / numworkers;

    /************************* master code *******************************/
    if (taskid == MASTER)
    {
        if (NX % numworkers != 0)
        {
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(1);
        }
        float u[2][NX][NY];
        inidat(NX, NY, &u[0]);
        // inidat(NX, NY, u);
        prtdat(NX, NY, &u[0], "initial.dat");

        offset = 0;
        // Send to workers
        for (i = 1; i <= numworkers; i++)
        {
            dest = i;
            if (dest == 1)
                up = NONE;
            else
                up = dest - 1;
            if (dest == numworkers)
                down = NONE;
            else
                down = dest + 1;

            MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&up, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&down, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&u[0][offset], rows * NY, MPI_FLOAT, dest, BEGIN,
                     MPI_COMM_WORLD);
            offset = offset + rows;
        }
        // Receive from works
        for (i = 1; i <= numworkers; i++)
        {
            source = i;
            msgtype = DONE;
            // Receive address of memory instead of data
            MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&u[0][offset], rows * NY, MPI_FLOAT, source,
                     msgtype, MPI_COMM_WORLD, &status);
        }
        prtdat(NX, NY, &u[0], "final.dat");
        MPI_Finalize();
    }

    /************************* workers code **********************************/
    if (taskid != MASTER)
    {
        float u[2][rows + 2][NY];
        // Initialize with zero
        for (iz = 0; iz < 2; iz++)
            for (ix = 0; ix < rows + 2; ix++)
                for (iy = 0; iy < NY; iy++)
                    u[iz][ix][iy] = 0.0;

        // Receive from master
        source = MASTER;
        msgtype = BEGIN;
        MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&u[0][1], rows * NY, MPI_FLOAT, source, msgtype,
                 MPI_COMM_WORLD, &status);

        iz = 0;
        for (it = 1; it <= STEPS; it++)
        {
            // Communicate with neighor works before computation
            if (up != NONE)
            {
                MPI_Send(&u[iz][1], NY, MPI_FLOAT, up,
                         DTAG, MPI_COMM_WORLD);
                source = up;
                msgtype = UTAG;
                MPI_Recv(&u[iz][0], NY, MPI_FLOAT, source,
                         msgtype, MPI_COMM_WORLD, &status);
            }
            if (down != NONE)
            {
                MPI_Send(&u[iz][rows], NY, MPI_FLOAT, down,
                         UTAG, MPI_COMM_WORLD);
                source = down;
                msgtype = DTAG;
                MPI_Recv(&u[iz][rows + 1], NY, MPI_FLOAT, source, msgtype,
                         MPI_COMM_WORLD, &status);
            }

            // Compute after sending and receiving data
            start = 1;
            end = rows;
            if (up == NONE)
            {
                start = 2;
            }
            if (down == NONE)
            {
                end--;
            }
            update(start, end, NY, &u[iz], &u[1 - iz]);
            iz = 1 - iz;
        }

        // Send final result to master
        MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&u[iz][1], rows * NY, MPI_FLOAT, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Finalize();
    }
}

void update(int start, int end, int ny, float *u1, float *u2)
{
    int ix, iy;
    for (ix = start; ix <= end; ix++)
        for (iy = 1; iy <= ny - 2; iy++)
            *(u2 + ix * ny + iy) = *(u1 + ix * ny + iy) +
                                   parms.cx * (*(u1 + (ix + 1) * ny + iy) +
                                               *(u1 + (ix - 1) * ny + iy) -
                                               2.0 * *(u1 + ix * ny + iy)) +
                                   parms.cy * (*(u1 + ix * ny + iy + 1) +
                                               *(u1 + ix * ny + iy - 1) -
                                               2.0 * *(u1 + ix * ny + iy));
}

void inidat(int nx, int ny, float *u)
{
    int ix, iy;

    for (ix = 0; ix <= nx - 1; ix++)
        for (iy = 0; iy <= ny - 1; iy++)
            if ((ix - nx / 2) * (ix - nx / 2) + (iy - ny / 2) * (iy - ny / 2) < 16)
            {
                *(u + ix * ny + iy) = 100.0;
            }
            else
            {
                *(u + ix * ny + iy) = 0.0;
            }
}

void prtdat(int nx, int ny, float *u1, char *fnam)
{
    int ix, iy;
    FILE *fp;

    fp = fopen(fnam, "w");
    for (iy = ny - 1; iy >= 0; iy--)
    {
        for (ix = 0; ix <= nx - 1; ix++)
        {
            fprintf(fp, "%8.1f", *(u1 + ix * ny + iy));
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}
