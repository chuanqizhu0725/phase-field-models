
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define N 2
#define NDX 128
#define NDY 128
#define NTH 8
#define PI 3.14159

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int rows;

int nstep = 4801;
int pstep = 1600;

double dx = 1.0;
double dtime = 0.02;
double gamma0 = 0.5;
double mobi = 1.0;
double delta = 5.0;

double A0, W0, M0, S0;

double Dl = 5.0;
double Ds = 0.01;

double DT = 0.0;
double temp0 = 2.0;
double gradT = 0.01;
double cl = 0.2;
// double cl = 0.5;

double sum0, lcount;
int fcount = 0;
double c0, cm0;
int i, j, im, ip, jm, jp, k, l;
int x00, y00;
double r0, r;

double calC01e(), calC1e(), calC02e(), calC2e();
double calDF10(), calDF20();

int main(void)
{
    double csta = dtime * Dl / dx / dx;
    double psta = dtime / dx / dx * mobi * A0;
    printf("----------------------------------------------\n");
    printf("Computation Started!");
    printf("concenration field stablity number is: %e\n", csta);
    printf("phase field stability number is: %e\n", psta);

    double(*mij)[N][N] = malloc(sizeof(*mij));
    double(*wij)[N][N] = malloc(sizeof(*wij));
    double(*aij)[N][N] = malloc(sizeof(*aij));
    double(*sij)[N][N] = malloc(sizeof(*sij));
    double(*phi)[N][NDX][NDY] = malloc(sizeof(*phi));
    double(*phi2)[N][NDX][NDY] = malloc(sizeof(*phi2));
    double(*cont)[NDX][NDY] = malloc(sizeof(*cont));
    double(*cont2)[NDX][NDY] = malloc(sizeof(*cont2));
    double(*conp)[N][NDX][NDY] = malloc(sizeof(*conp));
    double(*temp)[NDX][NDY] = malloc(sizeof(*temp));
    int(*phiNum)[NDX][NDY] = malloc(sizeof(*phiNum));
    int(*phiIdx)[N + 1][NDX][NDY] = malloc(sizeof(*phiIdx));

    A0 = 8.0 * delta * gamma0 / PI / PI;
    W0 = 4.0 * gamma0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);
    S0 = 0.5;

    // ------------------ Interface Initialization----------------------
    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            (*wij)[i][j] = W0;
            (*aij)[i][j] = A0;
            (*mij)[i][j] = M0;
            (*sij)[i][j] = S0;
            if (i == j)
            {
                (*wij)[i][j] = 0.0;
                (*aij)[i][j] = 0.0;
                (*mij)[i][j] = 0.0;
                (*sij)[i][j] = 0.0;
            }
        }
    }

    // ------------------ Temperature Initialization----------------------
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            (*temp)[i][j] = temp0 + i * gradT;
        }
    }

    // ------------------ Seeds Initialization   ----------------------
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            (*phi)[0][i][j] = 1.0;
            (*conp)[0][i][j] = cl;
            for (k = 1; k <= nm; k++)
            {

                (*phi)[k][i][j] = 0.0;
                (*conp)[k][i][j] = (k % 2 == 1) ? calC1e((*temp)[i][j]) : calC2e((*temp)[i][j]);
            }
        }
    }
    r0 = 15.0;
    for (k = 1; k <= nm; k++)
    {
        x00 = rand() % NDX;
        y00 = rand() % NDY;
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                r = sqrt((double)(((i - x00)) * (i - x00) + (j - y00) * (j - y00)));
                if (r <= r0)
                {
                    if (k % 2 == 1)
                    {
                        (*phi)[k][i][j] = 1.0;
                        (*conp)[k][i][j] = calC1e((*temp)[i][j]);
                        (*phi)[0][i][j] = 0.0;
                        (*conp)[0][i][j] = calC01e((*temp)[i][j]);
                        for (l = 0; l <= nm; l++)
                        {
                            if (l != k)
                            {
                                (*phi)[l][i][j] = 0.0;
                                (*conp)[l][i][j] = (l % 2 == 1) ? calC1e((*temp)[i][j]) : calC2e((*temp)[i][j]);
                            }
                        }
                    }
                    else if (k % 2 == 0)
                    {
                        (*phi)[k][i][j] = 1.0;
                        (*conp)[k][i][j] = calC2e((*temp)[i][j]);
                        (*phi)[0][i][j] = 0.0;
                        (*conp)[0][i][j] = calC02e((*temp)[i][j]);
                        for (l = 0; l <= nm; l++)
                        {
                            if (l != k)
                            {
                                (*phi)[l][i][j] = 0.0;
                                (*conp)[l][i][j] = (l % 2 == 1) ? calC1e((*temp)[i][j]) : calC2e((*temp)[i][j]);
                            }
                        }
                    }
                }
            }
        }
    }

    sum0 = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= nm; k++)
            {
                (*cont)[i][j] += (*phi)[k][i][j] * (*conp)[k][i][j];
            }
            sum0 += (*cont)[i][j];
        }
    }
    c0 = sum0 / NDX / NDY;
    cm0 = sum0;

    rows = NDX / NTH;

#pragma omp parallel num_threads(NTH)
    {

        int th_id = omp_get_thread_num();
        int offset, start, end, istep;
        int ix, iy, ixp, ixm, iyp, iym;
        int ii, jj, kk;
        int n1, n2, n3;
        double c1e, c01e, c2e, c02e;
        double dF, pddtt, sum1, sums, suml;
        double termiikk, termjjkk;
        double phidxx, phidyy;
        double cddtt;
        int phinum;
        double psum;

        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;
        istep = 0;

    //********************************** Computation start ******************************************
    start:;

        if ((istep % pstep == 0) && (th_id == 0))
        {
            FILE *stream; //ストリームのポインタ設定
            char buffer[30];
            sprintf(buffer, "data/phi/2d%d.csv", istep);
            stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    fprintf(stream, "%e   ", (*phi)[1][ix][iy]);
                    fprintf(stream, "\n");
                }
            }
            fclose(stream); //ファイルをクローズ

            FILE *streamc;
            char bufferc[30];
            sprintf(bufferc, "data/con/2d%d.csv", istep);
            streamc = fopen(bufferc, "a");

            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    fprintf(streamc, "%e   ", (*cont)[ix][iy]);
                    fprintf(streamc, "\n");
                }
            }
            fclose(streamc);

            FILE *streamt; //ストリームのポインタ設定
            char buffert[30];
            sprintf(buffert, "data/temp/2d%d.csv", istep);
            streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

            for (int ix = 0; ix <= ndmx; ix++)
            {
                fprintf(streamt, "%e   ", (*temp)[ix][0]);
                fprintf(streamt, "\n");
            }
            fclose(streamt); //ファイルをクローズ

            printf("%d steps have passed!\n", istep);
            printf("The nominal concnetration is %e\n", c0);
        }

        // Cooling rate
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                (*temp)[ix][iy] -= DT / 10000.0;
            }
        }

        // Collect phase field info
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (ix == ndmx)
                {
                    ixp = 0;
                }
                if (iy == ndmy)
                {
                    iyp = 0;
                }
                if (iy == 0)
                {
                    iym = ndmy;
                }

                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if (((*phi)[ii][ix][iy] > 0.0) ||
                        (((*phi)[ii][ix][iy] == 0.0) && ((*phi)[ii][ixp][iy] > 0.0) || ((*phi)[ii][ixm][iy] > 0.0) || ((*phi)[ii][ix][iyp] > 0.0) || ((*phi)[ii][ix][iym] > 0.0)))
                    {
                        phinum++;
                        (*phiIdx)[phinum][ix][iy] = ii;
                    }
                }
                (*phiNum)[ix][iy] = phinum;
            }
        }
#pragma omp barrier

        // Evolution Equation of Phase Fields
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (ix == ndmx)
                {
                    ixp = 0;
                }
                if (iy == ndmy)
                {
                    iyp = 0;
                }
                if (iy == 0)
                {
                    iym = ndmy;
                }
                for (n1 = 1; n1 <= (*phiNum)[ix][iy]; n1++)
                {
                    ii = (*phiIdx)[n1][ix][iy];
                    pddtt = 0.0;
                    for (n2 = 1; n2 <= (*phiNum)[ix][iy]; n2++)
                    {
                        jj = (*phiIdx)[n2][ix][iy];
                        sum1 = 0.0;
                        for (n3 = 1; n3 <= (*phiNum)[ix][iy]; n3++)
                        {
                            kk = (*phiIdx)[n3][ix][iy];

                            phidxx = ((*phi)[kk][ixp][iy] + (*phi)[kk][ixm][iy] - 2.0 * (*phi)[kk][ix][iy]) / (dx * dx);
                            phidyy = ((*phi)[kk][ix][iyp] + (*phi)[kk][ix][iym] - 2.0 * (*phi)[kk][ix][iy]) / (dx * dx);

                            termiikk = (*aij)[ii][kk] * phidxx;

                            termjjkk = (*aij)[jj][kk] * phidyy;

                            sum1 += 0.5 * (termiikk - termjjkk) + ((*wij)[ii][kk] - (*wij)[jj][kk]) * (*phi)[kk][ix][iy];
                        }
                        if ((ii == jj) || (ii > 0 && jj > 0))
                        {
                            dF = 0.0;
                        }
                        else if (ii % 2 == 1 && jj == 0)
                        {
                            dF = calDF10((*conp)[0][ix][iy], (*temp)[ix][iy], (*sij)[ii][jj]);
                        }
                        else if (ii == 0 && jj % 2 == 1)
                        {
                            dF = -calDF10((*conp)[0][ix][iy], (*temp)[ix][iy], (*sij)[ii][jj]);
                        }
                        else if (ii % 2 == 0 && jj == 0)
                        {
                            dF = calDF20((*conp)[0][ix][iy], (*temp)[ix][iy], (*sij)[ii][jj]);
                        }
                        else if (ii == 0 && jj % 2 == 0)
                        {
                            dF = -calDF20((*conp)[0][ix][iy], (*temp)[ix][iy], (*sij)[ii][jj]);
                        }
                        pddtt += -2.0 * (*mij)[ii][jj] / (double)((*phiNum)[ix][iy]) * (sum1 - 8.0 / PI * dF * sqrt((*phi)[ii][ix][iy] * (*phi)[jj][ix][iy]));
                    }
                    (*phi2)[ii][ix][iy] = (*phi)[ii][ix][iy] + pddtt * dtime;
                    if ((*phi2)[ii][ix][iy] >= 1.0)
                    {
                        (*phi2)[ii][ix][iy] = 1.0;
                    }
                    if ((*phi2)[ii][ix][iy] <= 0.0)
                    {
                        (*phi2)[ii][ix][iy] = 0.0;
                    }
                }
            } // i
        }

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (kk = 0; kk <= nm; kk++)
                {
                    (*phi)[kk][ix][iy] = (*phi2)[kk][ix][iy];
                }
            }
        }
        // local sum should be 1
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                psum = 0.0;
                for (kk = 0; kk <= nm; kk++)
                {
                    psum += (*phi)[kk][ix][iy];
                }
                for (kk = 0; kk <= nm; kk++)
                {
                    (*phi)[kk][ix][iy] = (*phi)[kk][ix][iy] / psum;
                }
            }
        }

#pragma omp barrier

        // Calculate the concentration field in solid and liqud phase
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                c1e = calC1e((*temp)[ix][iy]);
                c2e = calC2e((*temp)[ix][iy]);
                c01e = calC01e((*temp)[ix][iy]);
                c02e = calC02e((*temp)[ix][iy]);
                // for (kk = 1; kk <= (*phiNum)[ix][iy]; kk++)
                // {
                //     n1 = (*phiIdx)[kk][ix][iy];
                //     if (n1 != 0)
                //     {
                //         (*conp)[n1][ix][iy] = (n1 % 2 == 1) ? c1e : c2e;
                //     }
                // }
                if ((*phi)[0][ix][iy] == 0.0)
                {
                    (*conp)[0][ix][iy] = 0.0;
                    for (kk = 1; kk <= (*phiNum)[ix][iy]; kk++)
                    {
                        n2 = (*phiIdx)[kk][ix][iy];
                        if (n2 != 0)
                        {
                            (*conp)[0][ix][iy] += (n2 % 2 == 1) ? c01e * (*phi)[n2][ix][iy] : c02e * (*phi)[n2][ix][iy];
                        }
                    }
                }
                else if ((*phi)[0][ix][iy] > 0.0 && (*phi)[0][ix][iy] < 1.0)
                {
                    for (kk = 1; kk <= (*phiNum)[ix][iy]; kk++)
                    {
                        n1 = (*phiIdx)[kk][ix][iy];
                        if (n1 != 0)
                        {
                            (*conp)[n1][ix][iy] = (n1 % 2 == 1) ? c1e : c2e;
                        }
                    }
                    sum1 = (*cont)[ix][iy];
                    for (kk = 1; kk <= (*phiNum)[ix][iy]; kk++)
                    {
                        n2 = (*phiIdx)[kk][ix][iy];
                        if (n2 != 0)
                        {
                            sum1 -= (*phi)[n2][ix][iy] * (*conp)[n2][ix][iy];
                        }
                    }
                    (*conp)[0][ix][iy] = sum1 / (*phi)[0][ix][iy];
                    if ((*conp)[0][ix][iy] > 1.0)
                    {
                        (*conp)[0][ix][iy] = 1.0;
                    }
                    if ((*conp)[0][ix][iy] < 0.0)
                    {
                        (*conp)[0][ix][iy] = 0.0;
                    }
                }
                else if ((*phi)[0][ix][iy] == 1.0)
                {
                    (*conp)[0][ix][iy] = (*cont)[ix][iy];
                }
                (*cont)[ix][iy] = 0.0;
                for (kk = 1; kk <= (*phiNum)[ix][iy]; kk++)
                {
                    n3 = (*phiIdx)[kk][ix][iy];
                    (*cont)[ix][iy] += (*conp)[n3][ix][iy] * (*phi)[n3][ix][iy];
                }
                if ((*cont)[ix][iy] > 1.0)
                {
                    (*cont)[ix][iy] = 1.0;
                }
                if ((*cont)[ix][iy] < 0.0)
                {
                    (*cont)[ix][iy] = 0.0;
                }
            }
        }
#pragma omp barrier

        // Evolution Equation of Concentration field
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (ix == ndmx)
                {
                    ixp = 0;
                }
                if (iy == ndmy)
                {
                    iyp = 0;
                }
                if (iy == 0)
                {
                    iym = ndmy;
                }

                sums = 0.0;
                for (kk = 1; kk <= (*phiNum)[ix][iy]; kk++)
                {
                    n1 = (*phiIdx)[kk][ix][iy];
                    if (n1 != 0)
                    {
                        sums += 0.25 * (((*phi)[n1][ixp][iy] - (*phi)[n1][ixm][iy]) * ((*conp)[n1][ixp][iy] - (*conp)[n1][ixm][iy]) + ((*phi)[n1][ix][iyp] - (*phi)[n1][ix][iym]) * ((*conp)[n1][ix][iyp] - (*conp)[n1][ix][iym])) / dx / dx +
                                (*phi)[n1][ix][iy] * ((*conp)[n1][ixp][iy] + (*conp)[n1][ixm][iy] + (*conp)[n1][ix][iyp] + (*conp)[n1][ix][iym] - 4.0 * (*conp)[n1][ix][iy]) / dx / dx;
                    }
                }

                suml = 0.25 * (((*phi)[0][ixp][iy] - (*phi)[0][ixm][iy]) * ((*conp)[0][ixp][iy] - (*conp)[0][ixm][iy]) + ((*phi)[0][ix][iyp] - (*phi)[0][ix][iym]) * ((*conp)[0][ix][iyp] - (*conp)[0][ix][iym])) / dx / dx +
                       (*phi)[0][ix][iy] * ((*conp)[0][ixp][iy] + (*conp)[0][ixm][iy] + (*conp)[0][ix][iyp] + (*conp)[0][ix][iym] - 4.0 * (*conp)[0][ix][iy]) / dx / dx;

                cddtt = Ds * sums + Dl * suml;
                (*cont2)[ix][iy] = (*cont)[ix][iy] + cddtt * dtime;
            }
        }

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                (*cont)[ix][iy] = (*cont2)[ix][iy];
            }
        }
        fcount += 1;
#pragma omp barrier
        if (th_id == 0 && fcount == NTH)
        {
            // correction for mass conservation
            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    sum0 += (*cont)[ix][iy];
                    lcount += (*phi)[0][ix][iy];
                }
            }

            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    if ((*phi)[0][ix][iy] > 0.0)
                    {
                        (*conp)[0][ix][iy] = (*conp)[0][ix][iy] - (sum0 - cm0) / lcount * (*phi)[0][ix][iy];
                        (*cont)[ix][iy] = (*cont)[ix][iy] - (sum0 - cm0) / lcount * (*phi)[0][ix][iy];
                        if ((*conp)[0][ix][iy] > 1.0)
                        {
                            (*conp)[0][ix][iy] = 1.0;
                        }
                        if ((*conp)[0][ix][iy] < 0.0)
                        {
                            (*conp)[0][ix][iy] = 0.0;
                        }
                        if ((*cont)[ix][iy] > 1.0)
                        {
                            (*cont)[ix][iy] = 1.0;
                        }
                        if ((*cont)[ix][iy] < 0.0)
                        {
                            (*cont)[ix][iy] = 0.0;
                        }
                    }
                }
            }
            c0 = sum0 / NDX / NDY;
            sum0 = 0.0;
            lcount = 0.0;
            fcount = 0;

            // moving frame
        }
#pragma omp barrier

        istep = istep + 1;
        if (istep < nstep)
        {
            goto start;
        }

    end:;
    } // end of thread
    return 0;
}

double calC01e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml1 = -10.0;
    double kap1 = 0.2;
    double c01e;
    c01e = ce + (temp0 - Te) / ml1;
    return c01e;
}

double calC1e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml1 = -10.0;
    double kap1 = 0.2;
    double c1e;
    c1e = (ce + (temp0 - Te) / ml1) * kap1;
    return c1e;
}

double calC02e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml2 = 10.0;
    double kap2 = 0.2;
    double c02e;
    c02e = ce + (temp0 - Te) / ml2;
    return c02e;
}

double calC2e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml2 = 10.0;
    double kap2 = 0.2;
    double c2e;
    c2e = 1.0 - (1.0 - (ce + (temp0 - Te) / ml2)) * kap2;
    return c2e;
}

double calDF10(double con0, double temp0, double dS)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml1 = -10.0;
    double DF = ((con0 - ce) * ml1 + Te - temp0) * dS;
    return DF;
}

double calDF20(double con0, double temp0, double dS)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml2 = 10.0;
    double DF = ((con0 - ce) * ml2 + Te - temp0) * dS;
    return DF;
}