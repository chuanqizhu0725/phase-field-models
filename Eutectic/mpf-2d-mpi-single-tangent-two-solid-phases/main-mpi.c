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
    void inidat(),
        initialize(),
        compute(),
        savedata(),
        saveint(),
        update(),
        prtdat();

    int nm = N - 1;

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
    int aniso[N][N];

    int i, j, im, ip, jm, jp, k;
    int ii, jj, kk;
    int n1, n2, n3;
    int istep;

    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            aniso[i][j] = 0;
            thij[i][j] = 0.0;
            if ((i == 0) || (j == 0))
            {
                aniso[i][j] = 1;
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                aniso[i][j] = 0;
            }
        }
    }

    double c0, dc0;

    double phis, phil;

    double cddtt, dev1_s, dev2_s, dev1_l, dev2_l;

    double phis_ij, phis_ipj, phis_imj, phis_ijp, phis_ijm;
    double cons_ij, cons_ipj, cons_imj, cons_ijp, cons_ijm;

    double theta, theta0;
    double epsilon0;
    double termiikk, termjjkk;

    double phidx, phidy, phidxx, phidyy, phidxy;
    double ep, ep1p, ep2p;

    double dF;

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
    rows = NX / numworkers;

    /************************* master code *******************************/
    if (taskid == MASTER)
    {
        if (NX % numworkers != 0)
        {
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(1);
        }
        double
            con[2][NX][NY],
            con1[NX][NY],
            con2[NX][NY],
            con0[NX][NY];

        double phi[2][NX][NY][N];

        int phiNum[NY][NY];

        initialize(NX, NY, &phi[0], con0, con1, con2, con);
        savedata(NX, NY, &phi[0], 2, "initial.dat");

        float u[2][NX][NY];
        inidat(NX, NY, u);

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
            //// send phase fields
            MPI_Send(&phi[0][offset], rows * NY * N, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            //// send concentration fields
            MPI_Send(&con0[offset], rows * NY, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&con1[offset], rows * NY, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&con2[offset], rows * NY, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&con[0][offset], rows * NY, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);

            offset = offset + rows;
        }
        // Receive from workers
        for (i = 1; i <= numworkers; i++)
        {
            source = i;
            msgtype = DONE;
            MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&u[0][offset], rows * NY, MPI_FLOAT, source,
                     msgtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&phiNum[offset], rows * NY, MPI_INT, source,
                     msgtype, MPI_COMM_WORLD, &status);
            //// receive phase fields
            MPI_Recv(&phi[1][offset], rows * NY * N, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
            //// receive con fields
            MPI_Recv(&con0[offset], rows * NY, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&con1[offset], rows * NY, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&con2[offset], rows * NY, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&con[1][offset], rows * NY, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
        }
        prtdat(NX, NY, u, "final.dat");
        // savedata(NX, NY, &phi[1], "phi.dat");
        // saveint(NX, NY, phiNum, "phinum.dat");

        MPI_Finalize();
    }

    /************************* workers code **********************************/
    if (taskid != MASTER)
    {
        float u[2][rows + 2][NY];
        double
            con[2][rows + 2][NY],
            con1[rows + 2][NY],
            con2[rows + 2][NY],
            con0[rows + 2][NY];

        double phi[2][rows + 2][NY][N];

        int phiNum[rows + 2][NY];
        int phiIdx[N + 1][rows + 2][NY];
        // Initialize with zero
        for (ix = 0; ix < rows + 2; ix++)
        {
            for (iy = 0; iy < NY; iy++)
            {
                for (iz = 0; iz < 2; iz++)
                {
                    u[iz][ix][iy] = 0.0;
                    phi[iz][ix][iy][0] = 0.0;
                    phi[iz][ix][iy][1] = 0.0;
                    phi[iz][ix][iy][2] = 0.0;
                    con[iz][ix][iy] = 0.0;
                }
                con0[ix][iy] = 0.0;
                con1[ix][iy] = 0.0;
                con2[ix][iy] = 0.0;
            }
        }

        // Receive from master
        source = MASTER;
        msgtype = BEGIN;
        MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);

        MPI_Recv(&u[0][1], rows * NY, MPI_FLOAT, source, msgtype,
                 MPI_COMM_WORLD, &status);
        //// receive phase fields
        MPI_Recv(&phi[0][1], rows * NY * N, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        //// receive con fields
        MPI_Recv(&con0[1], rows * NY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&con1[1], rows * NY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&con2[1], rows * NY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&con[0][1], rows * NY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);

        iz = 0;
        for (it = 1; it <= STEPS; it++)
        {
            // Communicate with neighor works before computation
            if (up != NONE)
            {
                MPI_Send(&u[iz][1], NY, MPI_FLOAT, up,
                         DTAG, MPI_COMM_WORLD);
                //// send up boundaries of phase fields
                MPI_Send(&phi[iz][1], NY * N, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                //// send up boundaries of con fields
                MPI_Send(&con0[1], NY, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                MPI_Send(&con1[1], NY, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                MPI_Send(&con2[1], NY, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                MPI_Send(&con[iz][1], NY, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);

                source = up;
                msgtype = UTAG;
                MPI_Recv(&u[iz][0], NY, MPI_FLOAT, source,
                         msgtype, MPI_COMM_WORLD, &status);
                //// receive up boundaries of phase fields
                MPI_Recv(&phi[iz][0], NY * N, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                //// receive up boundaries of con fields
                MPI_Recv(&con0[0], NY, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&con1[0], NY, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&con2[0], NY, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&con[iz][0], NY, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
            }
            if (down != NONE)
            {
                MPI_Send(&u[iz][rows], NY, MPI_FLOAT, down,
                         UTAG, MPI_COMM_WORLD);
                //// send down boundaries of phase fields
                MPI_Send(&phi[iz][rows], NY * N, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                //// send down boundaries of con fields
                MPI_Send(&con0[rows], NY, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                MPI_Send(&con1[rows], NY, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                MPI_Send(&con2[rows], NY, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                MPI_Send(&con[iz][rows], NY, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);

                source = down;
                msgtype = DTAG;
                MPI_Recv(&u[iz][rows + 1], NY, MPI_FLOAT, source, msgtype,
                         MPI_COMM_WORLD, &status);
                //// receive down boundaries of phase fields
                MPI_Recv(&phi[iz][rows + 1], NY * N, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                //// receive down boundaries of con fields
                MPI_Recv(&con0[rows + 1], NY, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&con1[rows + 1], NY, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&con2[rows + 1], NY, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&con[iz][rows + 1], NY, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
            }

            // Compute after sending and receiving data
            start = 1;
            end = rows;
            if (up == NONE)
                start = 2;
            if (down == NONE)
                end--;
            update(start, end, NY, &u[iz], &u[1 - iz]);

            // compute(start, end, NX, NY,
            //         nm, phiNum, phiIdx,
            //         &phi[iz], &phi[1 - iz],
            //         aij, wij, mij,
            //         dx, dtime);

            iz = 1 - iz;
        }

        // Send final result to master
        MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&u[iz][1], rows * NY, MPI_FLOAT, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Send(&phiNum[1], rows * NY, MPI_INT, MASTER, DONE,
                 MPI_COMM_WORLD);
        //// send phase fields
        MPI_Send(&phi[iz][1], rows * NY * N, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        //// send con fields
        MPI_Send(&con0[1], rows * NY, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Send(&con1[1], rows * NY, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Send(&con2[1], rows * NY, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Send(&con[iz][1], rows * NY, MPI_DOUBLE, MASTER, DONE,
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

void initialize(int nx, int ny,
                double *phi,
                double *con0, double *con1, double *con2, double *con)
{
    int ix, iy;
    double sum1 = 0.0;

    for (ix = 0; ix <= nx - 1; ix++)
    {
        for (iy = 0; iy <= ny - 1; iy++)
        {
            if ((ix - nx / 2) * (ix - nx / 2) + (iy - ny / 2) * (iy - ny / 2) < 400)
            {
                *(phi + nx * ny + ix * ny + iy) = 1.0;
                *(con1 + ix * ny + iy) = 0.1;
                *(phi + nx * ny * 2 + ix * ny + iy) = 0.0;
                *(con2 + ix * ny + iy) = 0.8;
                *(phi + ix * ny + iy) = 0.0;
                *(con0 + ix * ny + iy) = 0.3;
            }
            else
            {
                *(phi + nx * ny + ix * ny + iy + 1) = 0.0;
                *(con1 + ix * ny + iy) = 0.1;
                *(phi + nx * ny * 2 + ix * ny + iy) = 0.0;
                *(con2 + ix * ny + iy) = 0.8;
                *(phi + ix * ny + iy) = 1.0;
                *(con0 + ix * ny + iy) = 0.2;
            }
            *(con + ix * ny + iy) = *(phi + ix * ny + iy + 1) * *(con1 + ix * ny + iy) + *(phi + ix * ny + iy + 2) * *(con2 + ix * ny + iy) + *(phi + ix * ny + iy) * *(con0 + ix * ny + iy);
        }
    }
}

void savedata(int nx, int ny, double *data, int p, char *fnam)
{
    int ix, iy;
    FILE *fp;

    fp = fopen(fnam, "w");
    for (ix = 0; ix <= nx - 1; ix++)
    {
        for (iy = 0; iy <= ny - 1; iy++)
        {
            fprintf(fp, "%e", *(data + ix * ny + iy + p * nx * ny));
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

// void compute(int start, int end, int nx, int ny,
//              int nm, int *phiNum, int *phiIdx,
//              double *phio, double *phin,
//              double *aij, double *wij, double *mij,
//              double dx, double dtime)
// {
//     // Collect local info of phase fields
//     int ix, iy,
//         phinum, ii, jj, kk;

//     for (ix = start; ix <= end; ix++)
//     {
//         for (iy = 1; iy <= ny - 2; iy++)
//         {
//             phinum = 0;
//             for (ii = 0; ii <= nm; ii++)
//             {
//                 if ((*(phio + ii * nx * ny + ix * ny + iy) > 0.0) ||
//                     ((*(phio + ii * nx * ny + ix * ny + iy) == 0.0) && (*(phio + ii * nx * ny + (ix + 1) * ny + iy) > 0.0) ||
//                      (*(phio + ii * nx * ny + (ix - 1) * ny + iy) > 0.0) ||
//                      (*(phio + ii * nx * ny + ix * ny + iy + 1) > 0.0) ||
//                      (*(phio + ii * nx * ny + ix * ny + iy - 1) > 0.0)))
//                 {
//                     phinum++;
//                     *(phiIdx + phinum * nx * ny + ix * ny + iy) = ii;
//                 }
//             }
//             *(phiNum + ix * ny + iy) = phinum;
//         }
//     }

// int n1, n2, n3;
// double dpdt, sum, termiikk, termjjkk;
// double dF;

// for (ix = start; ix <= end; ix++)
// {
//     for (iy = 1; iy <= ny - 2; iy++)
//     {

//         for (n1 = 1; n1 <= *(phiNum + ix * ny + iy); n1++)
//         {
//             ii = *(phiIdx + n1 * nx * ny + ix * ny + iy);
//             dpdt = 0.0;
//             for (n2 = 1; n2 <= *(phiNum + ix * ny + iy); n2++)
//             {
//                 jj = *(phiIdx + n2 * nx * ny + ix * ny + iy);
//                 sum = 0.0;
//                 for (n3 = 1; n3 <= *(phiNum + ix * ny + iy); n3++)
//                 {
//                     kk = *(phiIdx + n3 * nx * ny + ix * ny + iy);

//                     termiikk = *(aij + ii * (nm + 1) + kk) * (*(phio + kk * nx * ny + (ix + 1) * ny + iy) + *(phio + kk * nx * ny + (ix - 1) * ny + iy) + *(phio + kk * nx * ny + ix * ny + iy + 1) + *(phio + kk * nx * ny + ix * ny + iy - 1) - 4.0 * *(phio + kk * nx * ny + ix * ny + iy)) / (dx * dx);

//                     termjjkk = *(aij + jj * (nm + 1) + kk) * (*(phio + kk * nx * ny + (ix + 1) * ny + iy) + *(phio + kk * nx * ny + (ix - 1) * ny + iy) + *(phio + kk * nx * ny + ix * ny + iy + 1) + *(phio + kk * nx * ny + ix * ny + iy - 1) - 4.0 * *(phio + kk * nx * ny + ix * ny + iy)) / (dx * dx);

//                     sum += 0.5 * (termiikk - termjjkk) + (*(wij + ii * (nm + 1) + kk) - *(wij + jj * (nm + 1) + kk)) * *(phio + kk * nx * ny + ix * ny + iy);
//                 }
//                 if (ii == 1 && jj == 0)
//                 {
//                     dF = 2.0e6;
//                 }
//                 else if (ii == 0 && jj == 1)
//                 {
//                     dF = -2.0e6;
//                 }
//                 else
//                 {
//                     dF = 0.0;
//                 }
//                 dpdt += -2.0 * (*(mij + ii * (nm + 1) + jj) / (double)*(phiNum + ix * ny + iy)) * (sum - 8.0 / PI * dF * sqrt(*(phio + ii * nx * ny + ix * ny + iy) * *(phio + jj * nx * ny + ix * ny + iy)));
//             }
//             *(phin + ii * nx * ny + ix * ny + iy) = *(phio + ii * nx * ny + ix * ny + iy) + dpdt * dtime;
//             if (*(phin + ii * nx * ny + ix * ny + iy) >= 1.0)
//             {
//                 *(phin + ii * nx * ny + ix * ny + iy) = 1.0;
//             }
//             if (*(phin + ii * nx * ny + ix * ny + iy) <= 0.0)
//             {
//                 *(phin + ii * nx * ny + ix * ny + iy) = 0.0;
//             }
//         }
//     }
// }
// }

void saveint(int nx, int ny, int *data, char *fnam)
{
    int ix, iy;
    FILE *fp;

    fp = fopen(fnam, "w");
    for (ix = 0; ix <= nx - 1; ix++)
    {
        for (iy = 0; iy <= ny - 1; iy++)
        {
            fprintf(fp, "%d", *(data + ix * ny + iy));
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
