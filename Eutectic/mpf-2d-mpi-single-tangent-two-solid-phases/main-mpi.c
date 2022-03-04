#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 3
#define NX 100
#define NY 100
#define PI 3.14159
#define STEPS 100
#define BEGIN 1
#define UTAG 2
#define DTAG 3
#define NONE 0
#define DONE 4
#define MASTER 0
#define DRND(x) ((double)(x) / RAND_MAX * rand())

struct Parms
{
    float cx;
    float cy;
} parms = {0.1, 0.1};

int main(int argc, char *argv[])
{
    void inidat(), initialize(), savedata(), update(), prtdat();

    int nm = N - 1;
    int ndxm = NX - 1;

    int nstep = 50001;
    int pstep = 500;

    double dx = 1.0e-8;
    double dtime = 4.0e-12;

    double gamma0 = 0.1;
    double astre = 0.05;
    double delta = 5.0 * dx;

    double mobi = 4.2e-5;

    double A0 = 8.0 * delta * gamma0 / PI / PI;
    double W0 = 4.0 * gamma0 / delta;
    double M0 = mobi * PI * PI / (8.0 * delta);
    double F10 = 9.0e6;
    double F20 = 2.0e6;

    double Ds = 0.1e-6;
    double Dl = 0.2e-5;

    double temp = 600.0; // K
    double Te = 800.0;
    double ce = 0.4;
    double cl = 0.4;

    double ml1 = -2000.0;
    double kap1 = 0.25;
    double c01e = ce + (temp - Te) / ml1;
    double c1e = c01e * kap1;

    double ml2 = 2000.0;
    double kap2 = 2.5;
    double c02e = ce + (temp - Te) / ml2;
    double c2e = c02e * kap2;

    double mij[N][N], aij[N][N], wij[N][N], fij[N][N], thij[N][N];
    // bool anij[N][N];

    // int i, j;

    // for (i = 0; i <= nm; i++)
    // {
    //     for (j = 0; j <= nm; j++)
    //     {
    //         wij[i][j] = W0;
    //         aij[i][j] = A0;
    //         mij[i][j] = M0;
    //         anij[i][j] = false;
    //         thij[i][j] = 0.0;
    //         if ((i == 0) || (j == 0))
    //         {
    //             anij[i][j] = true;
    //         }
    //         if (i == j)
    //         {
    //             wij[i][j] = 0.0;
    //             aij[i][j] = 0.0;
    //             mij[i][j] = 0.0;
    //             anij[i][j] = false;
    //         }
    //     }
    // }

    int phinum;
    int phiNum[NX][NY];
    int phiIdx[N + 1][NX][NY];

    double c0, dc0;

    double phis, phil;

    double cddtt, dev1_s, dev2_s, dev1_l, dev2_l;

    double phis_ij, phis_ipj, phis_imj, phis_ijp, phis_ijm;
    double cons_ij, cons_ipj, cons_imj, cons_ijp, cons_ijm;

    int i, j, im, ip, jm, jp, k;
    int ii, jj, kk;
    int n1, n2, n3;
    int istep;

    double pddtt, sum1;

    double theta, theta0;
    double epsilon0;
    double termiikk, termjjkk;

    double phidx, phidy, phidxx, phidyy, phidxy;
    double ep, ep1p, ep2p;

    double dF;

    double con[2][NX][NY], con1[NX][NY], con2[NX][NY], con0[NX][NY];

    double phi[2][N][NX][NY];

    float u[2][NX][NY];
    int taskid,
        numworkers,
        numtasks,
        rows, offset,
        dest, source,
        up, down,
        msgtype,
        rc, start, end,
        ix, iy, iz, it;

    MPI_Status status;

    // Allocate taskid to each core (cores = tasks = master + workers)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    numworkers = numtasks - 1;

    /************************* master code *******************************/
    if (taskid == MASTER)
    {
        if (NX % numworkers != 0)
        {
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(1);
        }

        initialize(NX, NY, &phi[0][0], &phi[0][1], &phi[0][2], con0, con1, con2, con);
        inidat(NX, NY, &u[0]);
        savedata(NX, NY, &phi[0][1], "initial.dat");

        rows = NX / numworkers;
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
            MPI_Send(&u[0][offset][0], rows * NY, MPI_FLOAT, dest, BEGIN,
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
            MPI_Recv(&u[0][offset][0], rows * NY, MPI_FLOAT, source,
                     msgtype, MPI_COMM_WORLD, &status);
        }
        prtdat(NX, NY, &u[0][0][0], "final.dat");
        MPI_Finalize();
    }

    /************************* workers code **********************************/
    if (taskid != MASTER)
    {
        // Initialize with zero
        for (iz = 0; iz < 2; iz++)
            for (ix = 0; ix < NX; ix++)
                for (iy = 0; iy < NY; iy++)
                    u[iz][ix][iy] = 0.0;

        // Receive from master
        source = MASTER;
        msgtype = BEGIN;
        MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&u[0][offset][0], rows * NY, MPI_FLOAT, source, msgtype,
                 MPI_COMM_WORLD, &status);

        iz = 0;
        for (it = 1; it <= STEPS; it++)
        {
            // Communicate with neighor works before computation
            if (up != NONE)
            {
                MPI_Send(&u[iz][offset][0], NY, MPI_FLOAT, up,
                         DTAG, MPI_COMM_WORLD);
                source = up;
                msgtype = UTAG;
                MPI_Recv(&u[iz][offset - 1][0], NY, MPI_FLOAT, source,
                         msgtype, MPI_COMM_WORLD, &status);
            }
            if (down != NONE)
            {
                MPI_Send(&u[iz][offset + rows - 1][0], NY, MPI_FLOAT, down,
                         UTAG, MPI_COMM_WORLD);
                source = down;
                msgtype = DTAG;
                MPI_Recv(&u[iz][offset + rows][0], NY, MPI_FLOAT, source, msgtype,
                         MPI_COMM_WORLD, &status);
            }

            // Compute after sending and receiving data
            start = offset;
            end = offset + rows - 1;
            if (offset == 0)
                start = 1;
            if ((offset + rows) == NX)
                end--;
            update(start, end, NY, &u[iz][0][0], &u[1 - iz][0][0]);
            iz = 1 - iz;
        }

        // Send final result to master
        MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&u[iz][offset][0], rows * NY, MPI_FLOAT, MASTER, DONE,
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

void initialize(int nx, int ny, double *phi0, double *phi1, double *phi2, double *con0, double *con1, double *con2, double *con)
{
    int ix, iy;
    double sum1 = 0.0;

    for (ix = 0; ix <= nx - 1; ix++)
    {
        for (iy = 0; iy <= ny - 1; iy++)
        {
            if ((ix - nx / 2) * (ix - nx / 2) + (iy - ny / 2) * (iy - ny / 2) < 400)
            {
                *(phi1 + ix * ny + iy) = 1.0;
                *(con1 + ix * ny + iy) = 0.1;
                *(phi2 + ix * ny + iy) = 0.0;
                *(con2 + ix * ny + iy) = 0.8;
                *(phi0 + ix * ny + iy) = 0.0;
                *(con0 + ix * ny + iy) = 0.3;
            }
            else
            {
                *(phi1 + ix * ny + iy) = 0.0;
                *(con1 + ix * ny + iy) = 0.1;
                *(phi2 + ix * ny + iy) = 0.0;
                *(con2 + ix * ny + iy) = 0.8;
                *(phi0 + ix * ny + iy) = 1.0;
                *(con0 + ix * ny + iy) = 0.2;
            }
            *(con + ix * ny + iy) = *(phi1 + ix * ny + iy) * *(con1 + ix * ny + iy) + *(phi2 + ix * ny + iy) * *(con2 + ix * ny + iy) + *(phi0 + ix * ny + iy) * *(con0 + ix * ny + iy);
        }
    }
}

void savedata(int nx, int ny, double *data, char *fnam)
{
    int ix, iy;
    FILE *fp;

    fp = fopen(fnam, "w");
    for (ix = 0; ix <= nx - 1; ix++)
    {
        for (iy = 0; iy <= ny - 1; iy++)
        {
            fprintf(fp, "%e", *(data + ix * ny + iy));
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
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
            fprintf(fp, "%e", *(u1 + ix * ny + iy));
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
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