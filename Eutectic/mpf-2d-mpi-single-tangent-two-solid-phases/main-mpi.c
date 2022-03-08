#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define N 3
#define NX 50
#define NY 50
#define PI 3.14159
#define STEPS 1
#define BEGIN 1
#define UTAG 2
#define DTAG 3
#define NONE 0
#define DONE 4
#define MASTER 0
#define DRND(x) ((double)(x) / RAND_MAX * rand())

int main(int argc, char *argv[])
{
    void savedata(), saveint();

    int nm = N - 1;
    int ndmx = NX - 1;
    int ndmy = NY - 1;

    int i, j, k, l, ii, jj, kk, ll; //整数
    int ip, im, jp, jm;             //整数
    int n1, n2, n3;                 //整数

    int nstep; //計算カウント数の最大値（計算終了カウント）

    double M0;               //粒界の易動度
    double W0;               //ペナルティー項の係数
    double A0;               //勾配エネルギー係数
    double F0;               //粒界移動の駆動力
    double sum1, sum2, sum3; //各種の和の作業変数
    double pddtt;            //フェーズフィールドの時間変化率

    double gamma0;      //粒界エネルギ密度
    double RR = 8.3145; //ガス定数

    double L = 2000.0;
    double dx = L / (double)NX * 1.0e-9; //差分プロック１辺の長さ(m)
    double dtime = 5.0;
    double temp = 1000.0;
    double vm0 = 7.0e-6;
    double delta = dx * 7.0;
    double mobi = 1.0;
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
    W0 = 4.0 * gamma0 / delta;           //ペナルティー項の係数[式(4.40)]
    M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]
    F0 = 50.0 / RR / temp;               //粒界移動の駆動力

    double aij[N][N]; //勾配エネルギー係数
    double wij[N][N]; //ペナルティー項の係数
    double mij[N][N]; //粒界の易動度
    double fij[N][N]; //粒界移動の駆動力

    for (i = 1; i <= nm; i++)
    {
        for (j = 1; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            fij[i][j] = 0.0;
            if ((i == nm) || (j == nm))
            {
                fij[i][j] = F0;
            }
            if (i > j)
            {
                fij[i][j] = -fij[i][j];
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                fij[i][j] = 0.0;
            }
        }
    }

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

        double phi[2][N][NX][NY];
        double intphi[NX][NY];
        int phiNum[NX][NY];
        double phi0[NX][NY];

        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                if ((i - NX / 2) * (i - NX / 2) + (j - NY / 2) * (j - NY / 2) <= 400)
                {
                    phi[0][1][i][j] = 1.0;
                    phi[0][2][i][j] = 0.0;
                    phi[0][0][i][j] = 0.0;
                }
                else
                {
                    phi[0][1][i][j] = 0.0;
                    phi[0][2][i][j] = 0.0;
                    phi[0][0][i][j] = 1.0;
                }
            }
        }

        savedata(NX, NY, &phi[0][1], "initialphi.dat");

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
            //// send phase fields
            MPI_Send(&phi[0][0][offset], rows * NY, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&phi[0][1][offset], rows * NY, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&phi[0][2][offset], rows * NY, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);

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
            MPI_Recv(&phiNum[offset], rows * NY, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
            //// receive phase fields
            MPI_Recv(&intphi[offset], rows * NY, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&phi0[offset], rows * NY, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
        }
        savedata(NX, NY, intphi, "intphi.dat");
        savedata(NX, NY, phi0, "phi0.dat");
        saveint(NX, NY, phiNum, "phinum.dat");

        MPI_Finalize();
    }

    /************************* workers code **********************************/
    if (taskid != MASTER)
    {

        double phi[2][N][rows + 2][NY];
        int phiNum[rows + 2][NY];
        int phiIdx[N + 1][rows + 2][NY];
        double intphi[rows + 2][NY];
        double phi0[rows + 2][NY];
        // Initialize with zero
        for (ix = 0; ix <= rows + 1; ix++)
        {
            for (iy = 0; iy <= NY - 1; iy++)
            {
                for (iz = 0; iz <= 1; iz++)
                {
                    phi[iz][0][ix][iy] = 0.0;
                    phi[iz][1][ix][iy] = 0.0;
                    phi[iz][2][ix][iy] = 0.0;
                }
                phiNum[ix][iy] = 0;
            }
        }

        // Receive from master
        source = MASTER;
        msgtype = BEGIN;
        MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        //// receive phase fields
        MPI_Recv(&phi[0][0][1], rows * NY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&phi[0][1][1], rows * NY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&phi[0][2][1], rows * NY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);

        iz = 0;
        for (it = 1; it <= STEPS; it++)
        {
            // Communicate with neighor works before computation
            if (up != NONE)
            {
                //// send up boundaries of phase fields
                MPI_Send(&phi[iz][0][1], NY, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                MPI_Send(&phi[iz][1][1], NY, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                MPI_Send(&phi[iz][2][1], NY, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);

                source = up;
                msgtype = UTAG;
                //// receive up boundaries of phase fields
                MPI_Recv(&phi[iz][0][0], NY, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&phi[iz][1][0], NY, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&phi[iz][2][0], NY, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
            }
            if (down != NONE)
            {
                //// send down boundaries of phase fields
                MPI_Send(&phi[iz][0][rows], NY, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                MPI_Send(&phi[iz][1][rows], NY, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                MPI_Send(&phi[iz][2][rows], NY, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);

                source = down;
                msgtype = DTAG;
                //// receive down boundaries of phase fields
                MPI_Recv(&phi[iz][0][rows + 1], NY, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&phi[iz][1][rows + 1], NY, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&phi[iz][2][rows + 1], NY, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
            }

            // Compute after sending and receiving data
            start = 1;
            end = rows;
            // update(start, end, NY, &u[iz], &u[1 - iz]);

            for (i = start; i <= end; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    if (up == NONE && i == 1)
                    {
                        im = 1;
                    }
                    if (down == NONE && i == rows)
                    {
                        ip = rows;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }

                    int phinum = 0;
                    for (ii = 0; ii <= nm; ii++)
                    {
                        if ((phi[0][ii][i][j] > 0.0) ||
                            ((phi[0][ii][i][j] == 0.0) && (phi[0][ii][ip][j] > 0.0) ||
                             (phi[0][ii][im][j] > 0.0) ||
                             (phi[0][ii][i][jp] > 0.0) ||
                             (phi[0][ii][i][jm] > 0.0)))
                        {
                            phinum++;
                            phiIdx[phinum][i][j] = ii;
                        }
                    }
                    phiNum[i][j] = phinum;
                }
            }

            // Evolution Equations
            for (i = start; i <= end; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    if (up == NONE && i == 1)
                    {
                        im = 1;
                    }
                    if (down == NONE && i == rows)
                    {
                        ip = rows;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }

                    for (n1 = 1; n1 <= phiNum[i][j]; n1++)
                    {
                        ii = phiIdx[n1][i][j];
                        pddtt = 0.0;
                        for (n2 = 1; n2 <= phiNum[i][j]; n2++)
                        {
                            jj = phiIdx[n2][i][j];
                            sum1 = 0.0;
                            for (n3 = 1; n3 <= phiNum[i][j]; n3++)
                            {
                                kk = phiIdx[n3][i][j];
                                sum1 += 0.5 * (aij[ii][kk] - aij[jj][kk]) * (phi[iz][kk][ip][j] + phi[iz][kk][im][j] + phi[iz][kk][i][jp] + phi[iz][kk][i][jm] - 4.0 * phi[iz][kk][i][j]) + (wij[ii][kk] - wij[jj][kk]) * phi[iz][kk][i][j]; //[式(4.31)の一部]
                            }
                            pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j] * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt(phi[iz][ii][i][j] * phi[iz][jj][i][j]));
                            //フェーズフィールドの発展方程式[式(4.31)]
                        }
                        if (pddtt >= 10e-20)
                        {
                            printf("Computation is going!\n");
                            printf("%e\n", pddtt);
                        }
                        phi[1 - iz][ii][i][j] = phi[iz][ii][i][j] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
                        if (phi[1 - iz][ii][i][j] >= 1.0)
                        {
                            phi[1 - iz][ii][i][j] = 1.0;
                        } //フェーズフィールドの変域補正
                        if (phi[1 - iz][ii][i][j] <= 0.0)
                        {
                            phi[1 - iz][ii][i][j] = 0.0;
                        }
                    }
                } // j
            }     // i
            iz = 1 - iz;
        }

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (ii = 0; ii <= nm; ii++)
                {
                    intphi[i][j] += phi[iz][ii][i][j] * phi[iz][ii][i][j];
                }
                phi0[i][j] = phi[iz][0][i][j];
            }
        }
        // Send final result to master
        MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&phiNum[1], rows * NY, MPI_INT, MASTER, DONE,
                 MPI_COMM_WORLD);
        //// send phase fields
        MPI_Send(&intphi[1], rows * NY, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Send(&phi0[1], rows * NY, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Finalize();
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

// void prtdat(int nx, int ny, float *u1, char *fnam)
// {
//     int ix, iy;
//     FILE *fp;

//     fp = fopen(fnam, "w");
//     for (iy = ny - 1; iy >= 0; iy--)
//     {
//         for (ix = 0; ix <= nx - 1; ix++)
//         {
//             fprintf(fp, "%e", *(u1 + ix * ny + iy));
//             fprintf(fp, "\n");
//         }
//     }
//     fclose(fp);
// }

// void inidat(int nx, int ny, float *u)
// {
//     int ix, iy;

//     for (ix = 0; ix <= nx - 1; ix++)
//         for (iy = 0; iy <= ny - 1; iy++)
//             if ((ix - nx / 2) * (ix - nx / 2) + (iy - ny / 2) * (iy - ny / 2) < 16)
//             {
//                 *(u + ix * ny + iy) = 100.0;
//             }
//             else
//             {
//                 *(u + ix * ny + iy) = 0.0;
//             }
// }
