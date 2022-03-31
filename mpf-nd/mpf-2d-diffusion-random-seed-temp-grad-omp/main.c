
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define N 9
#define NDX 64
#define NDY 64
#define NTH 8
#define PI 3.14159

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;

int nstep = 801;
int pstep = 400;

double dx = 1.0;
double dtime = 0.02;
double gamma0 = 0.5;
double mobi = 1.0;
double delta = 5.0;

double A0, W0, M0, S0;

double Dl = 5.0;
double Ds = 0.01;

double DT = 0.0;
double temp0 = -1.0;
double gradT = 0.00;
// double cl = 0.2;
double cl = 0.5;

double c0, cm0;
double sum0, lcount;

int x00, y00;
double r0, r;

double calC01e(), calC1e(), calC02e(), calC2e();
double calDF10(), calDF20();

int main(void)
{
    A0 = 8.0 * delta * gamma0 / PI / PI;
    W0 = 4.0 * gamma0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);
    S0 = 0.5;

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

    // ------------------ Interface Initialization----------------------
    for (int ni = 0; ni <= nm; ni++)
    {
        for (int nj = 0; nj <= nm; nj++)
        {
            (*wij)[ni][nj] = W0;
            (*aij)[ni][nj] = A0;
            (*mij)[ni][nj] = M0;
            (*sij)[ni][nj] = S0;
            if (ni == nj)
            {
                (*wij)[ni][nj] = 0.0;
                (*aij)[ni][nj] = 0.0;
                (*mij)[ni][nj] = 0.0;
                (*sij)[ni][nj] = 0.0;
            }
        }
    }

    // ------------------ Temperature Initialization----------------------
    for (int ix = 0; ix <= ndmx; ix++)
    {
        for (int iy = 0; iy <= ndmy; iy++)
        {
            (*temp)[ix][iy] = temp0 + ix * gradT;
        }
    }

    // ------------------ Seeds Initialization   ----------------------
    for (int ix = 0; ix <= ndmx; ix++)
    {
        for (int iy = 0; iy <= ndmy; iy++)
        {
            (*phi)[0][ix][iy] = 1.0;
            (*conp)[0][ix][iy] = cl;
            for (int nk = 1; nk <= nm; nk++)
            {

                (*phi)[nk][ix][iy] = 0.0;
                (*conp)[nk][ix][iy] = (nk % 2 == 1) ? calC1e((*temp)[ix][iy]) : calC2e((*temp)[ix][iy]);
            }
        }
    }
    r0 = 5.0;
    for (int ni = 1; ni <= nm; ni++)
    {
        x00 = rand() % NDX;
        y00 = rand() % NDY;
        for (int ix = 0; ix <= ndmx; ix++)
        {
            for (int iy = 0; iy <= ndmy; iy++)
            {
                r = sqrt((double)(((ix - x00)) * (ix - x00) + (iy - y00) * (iy - y00)));
                if (r <= r0)
                {
                    if (ni % 2 == 1)
                    {
                        (*phi)[ni][ix][iy] = 1.0;
                        (*conp)[ni][ix][iy] = calC1e((*temp)[ix][iy]);
                        (*phi)[0][ix][iy] = 0.0;
                        (*conp)[0][ix][iy] = calC01e((*temp)[ix][iy]);
                        for (int nk = 0; nk <= nm; nk++)
                        {
                            if (nk != ni)
                            {
                                (*phi)[nk][ix][iy] = 0.0;
                                (*conp)[nk][ix][iy] = (nk % 2 == 1) ? calC1e((*temp)[ix][iy]) : calC2e((*temp)[ix][iy]);
                            }
                        }
                    }
                    else if (ni % 2 == 0)
                    {
                        (*phi)[ni][ix][iy] = 1.0;
                        (*conp)[ni][ix][iy] = calC2e((*temp)[ix][iy]);
                        (*phi)[0][ix][iy] = 0.0;
                        (*conp)[0][ix][iy] = calC02e((*temp)[ix][iy]);
                        for (int nk = 0; nk <= nm; nk++)
                        {
                            if (nk != ni)
                            {
                                (*phi)[nk][ix][iy] = 0.0;
                                (*conp)[nk][ix][iy] = (nk % 2 == 1) ? calC1e((*temp)[ix][iy]) : calC2e((*temp)[ix][iy]);
                            }
                        }
                    }
                }
            }
        }
    }

    sum0 = 0.0;
    for (int ix = 0; ix <= ndmx; ix++)
    {
        for (int iy = 0; iy <= ndmy; iy++)
        {
            for (int nk = 0; nk <= nm; nk++)
            {
                (*cont)[ix][iy] += (*phi)[nk][ix][iy] * (*conp)[nk][ix][iy];
            }
            sum0 += (*cont)[ix][iy];
        }
    }
    c0 = sum0 / NDX / NDY;
    cm0 = sum0;

    //********************************** Computation start ******************************************

    int rows = NDX / NTH;
    sum0 = 0.0;
    lcount = 0.0;

#pragma omp parallel num_threads(NTH)
    {
        int i, j, im, ip, jm, jp;
        int ii, jj, kk;
        int n1, n2, n3, phinum;
        double c1e, c01e, c2e, c02e;
        double dc0, dcm0;
        double cddtt;
        double phidxx, phidyy;
        double termiikk, termjjkk;
        double dF, pddtt, sum1, sums, suml;

        int th_id, offset, start, end, istep;

        th_id = omp_get_thread_num();
        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;
        istep = 0;

    start:;

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                if (i == ndmx)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndmx;
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }

                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if (((*phi)[ii][i][j] > 0.0) ||
                        (((*phi)[ii][i][j] == 0.0) && ((*phi)[ii][ip][j] > 0.0) || ((*phi)[ii][im][j] > 0.0) || ((*phi)[ii][i][jp] > 0.0) || ((*phi)[ii][i][jm] > 0.0)))
                    {
                        phinum++;
                        (*phiIdx)[phinum][i][j] = ii;
                    }
                }
                (*phiNum)[i][j] = phinum;
            }
        }

        // Calculate the concentration field in solid and liqud phase
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                c1e = calC1e((*temp)[i][j]);
                c2e = calC2e((*temp)[i][j]);
                c01e = calC01e((*temp)[i][j]);
                c02e = calC02e((*temp)[i][j]);
                for (kk = 1; kk <= (*phiNum)[i][j]; kk++)
                {
                    n1 = (*phiIdx)[kk][i][j];
                    if (n1 != 0)
                    {
                        (*conp)[n1][i][j] = (n1 % 2 == 1) ? c1e : c2e;
                    }
                }
                if ((*phi)[0][i][j] == 0.0)
                {
                    (*conp)[0][i][j] = 0.0;
                    for (kk = 1; kk <= (*phiNum)[i][j]; kk++)
                    {
                        n2 = (*phiIdx)[kk][i][j];
                        if (n2 != 0)
                        {
                            (*conp)[0][i][j] += (n2 % 2 == 1) ? c01e * (*phi)[n2][i][j] : c02e * (*phi)[n2][i][j];
                        }
                    }
                }
                else if ((*phi)[0][i][j] > 0.0 && (*phi)[0][i][j] < 1.0)
                {
                    sum1 = (*cont)[i][j];
                    for (kk = 1; kk <= (*phiNum)[i][j]; kk++)
                    {
                        n2 = (*phiIdx)[kk][i][j];
                        if (n2 != 0)
                        {
                            sum1 -= (*phi)[n2][i][j] * (*conp)[n2][i][j];
                        }
                    }
                    (*conp)[0][i][j] = sum1 / (*phi)[0][i][j];
                    if ((*conp)[0][i][j] > 1.0)
                    {
                        (*conp)[0][i][j] = 1.0;
                    }
                    if ((*conp)[0][i][j] < 0.0)
                    {
                        (*conp)[0][i][j] = 0.0;
                    }
                }
                else if ((*phi)[0][i][j] == 1.0)
                {
                    (*conp)[0][i][j] = (*cont)[i][j];
                }
                (*cont)[i][j] = 0.0;
                for (kk = 1; kk <= (*phiNum)[i][j]; kk++)
                {
                    n3 = (*phiIdx)[kk][i][j];
                    (*cont)[i][j] += (*conp)[n3][i][j] * (*phi)[n3][i][j];
                }
                if ((*cont)[i][j] > 1.0)
                {
                    (*cont)[i][j] = 1.0;
                }
                if ((*cont)[i][j] < 0.0)
                {
                    (*cont)[i][j] = 0.0;
                }
            }
        }

        // Evolution Equation of Concentration field
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                if (i == ndmx)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndmx;
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }

                sums = 0.0;
                for (kk = 1; kk <= (*phiNum)[i][j]; kk++)
                {
                    n1 = (*phiIdx)[kk][i][j];
                    if (n1 != 0)
                    {
                        sums += 0.25 * (((*phi)[n1][ip][j] - (*phi)[n1][im][j]) * ((*conp)[n1][ip][j] - (*conp)[n1][im][j]) + ((*phi)[n1][i][jp] - (*phi)[n1][i][jm]) * ((*conp)[n1][i][jp] - (*conp)[n1][i][jm])) / dx / dx +
                                (*phi)[n1][i][j] * ((*conp)[n1][ip][j] + (*conp)[n1][im][j] + (*conp)[n1][i][jp] + (*conp)[n1][i][jm] - 4.0 * (*conp)[n1][i][j]) / dx / dx;
                    }
                }

                suml = 0.25 * (((*phi)[0][ip][j] - (*phi)[0][im][j]) * ((*conp)[0][ip][j] - (*conp)[0][im][j]) + ((*phi)[0][i][jp] - (*phi)[0][i][jm]) * ((*conp)[0][i][jp] - (*conp)[0][i][jm])) / dx / dx +
                       (*phi)[0][i][j] * ((*conp)[0][ip][j] + (*conp)[0][im][j] + (*conp)[0][i][jp] + (*conp)[0][i][jm] - 4.0 * (*conp)[0][i][j]) / dx / dx;

                cddtt = Ds * sums + Dl * suml;
                (*cont2)[i][j] = (*cont)[i][j] + cddtt * dtime; //濃度場の時間発展(陽解法)
                                                                // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
            }
        }

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                (*cont)[i][j] = (*cont2)[i][j];
                sum0 += (*cont)[i][j];
            }
        }

        // #pragma omp barrier
        //         // calculate extra mass
        //         c0 = sum0 / NDX / NDY;
        //         dcm0 = sum0 - cm0;
        //         // collect mass in liquid
        //         for (i = start; i <= end; i++)
        //         {
        //             for (j = 0; j <= ndmy; j++)
        //             {
        //                 if ((*phi)[0][i][j] > 0.0)
        //                 {
        //                     lcount += (*phi)[0][i][j];
        //                 }
        //             }
        //         }
        // #pragma omp barrier
        //         // correction for mass conservation
        //         for (i = start; i <= end; i++)
        //         {
        //             for (j = 0; j <= ndmy; j++)
        //             {
        //                 if ((*phi)[0][i][j] > 0.0)
        //                 {
        //                     (*conp)[0][i][j] = (*conp)[0][i][j] - dcm0 / lcount * (*phi)[0][i][j];
        //                     (*cont)[i][j] = (*cont)[i][j] - dcm0 / lcount * (*phi)[0][i][j];
        //                     if ((*conp)[0][i][j] > 1.0)
        //                     {
        //                         (*conp)[0][i][j] = 1.0;
        //                     }
        //                     if ((*conp)[0][i][j] < 0.0)
        //                     {
        //                         (*conp)[0][i][j] = 0.0;
        //                     }
        //                 }
        //             }
        //         }

        // Evolution Equation of Phase Fields
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                if (i == ndmx)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndmx;
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                for (n1 = 1; n1 <= (*phiNum)[i][j]; n1++)
                {
                    ii = (*phiIdx)[n1][i][j];
                    pddtt = 0.0;
                    for (n2 = 1; n2 <= (*phiNum)[i][j]; n2++)
                    {
                        jj = (*phiIdx)[n2][i][j];
                        sum1 = 0.0;
                        for (n3 = 1; n3 <= (*phiNum)[i][j]; n3++)
                        {
                            kk = (*phiIdx)[n3][i][j];

                            phidxx = ((*phi)[kk][ip][j] + (*phi)[kk][im][j] - 2.0 * (*phi)[kk][i][j]) / (dx * dx);
                            phidyy = ((*phi)[kk][i][jp] + (*phi)[kk][i][jm] - 2.0 * (*phi)[kk][i][j]) / (dx * dx);

                            termiikk = (*aij)[ii][kk] * phidxx;

                            termjjkk = (*aij)[jj][kk] * phidyy;

                            sum1 += 0.5 * (termiikk - termjjkk) + ((*wij)[ii][kk] - (*wij)[jj][kk]) * (*phi)[kk][i][j];
                        }
                        if ((ii == jj) || (ii > 0 && jj > 0))
                        {
                            dF = 0.0;
                        }
                        else if (ii % 2 == 1 && jj == 0)
                        {
                            dF = calDF10((*conp)[0][i][j], (*temp)[i][j], (*sij)[ii][jj]);
                        }
                        else if (ii == 0 && jj % 2 == 1)
                        {
                            dF = -calDF10((*conp)[0][i][j], (*temp)[i][j], (*sij)[ii][jj]);
                        }
                        else if (ii % 2 == 0 && jj == 0)
                        {
                            dF = calDF20((*conp)[0][i][j], (*temp)[i][j], (*sij)[ii][jj]);
                        }
                        else if (ii == 0 && jj % 2 == 0)
                        {
                            dF = -calDF20((*conp)[0][i][j], (*temp)[i][j], (*sij)[ii][jj]);
                        }
                        pddtt += -2.0 * (*mij)[ii][jj] / (double)((*phiNum)[i][j]) * (sum1 - 8.0 / PI * dF * sqrt((*phi)[ii][i][j] * (*phi)[jj][i][j]));
                    }
                    (*phi2)[ii][i][j] = (*phi)[ii][i][j] + pddtt * dtime;
                    if ((*phi2)[ii][i][j] >= 1.0)
                    {
                        (*phi2)[ii][i][j] = 1.0;
                    }
                    if ((*phi2)[ii][i][j] <= 0.0)
                    {
                        (*phi2)[ii][i][j] = 0.0;
                    }
                }
            } // i
        }

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (kk = 0; kk <= nm; kk++)
                {
                    (*phi)[kk][i][j] = (*phi2)[kk][i][j];
                }
            }
        }

        // cooling rate
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                (*temp)[i][j] -= DT / 10000.0;
            }
        }

#pragma omp barrier

        if (th_id == 0)
        {
            sum0 = 0.0;
            lcount = 0.0;
        }

        if ((th_id == 0) && (istep % pstep == 0))
        {
            FILE *stream; //ストリームのポインタ設定
            char buffer[30];
            sprintf(buffer, "data/phi/2d%d.csv", istep);
            stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

            for (int i = 0; i <= ndmx; i++)
            {
                for (int j = 0; j <= ndmy; j++)
                {
                    fprintf(stream, "%e   ", (*phi)[1][i][j]);
                    fprintf(stream, "\n");
                }
            }
            fclose(stream); //ファイルをクローズ

            FILE *streamc;
            char bufferc[30];
            sprintf(bufferc, "data/con/2d%d.csv", istep);
            streamc = fopen(bufferc, "a");

            for (int i = 0; i <= ndmx; i++)
            {
                for (int j = 0; j <= ndmy; j++)
                {
                    fprintf(streamc, "%e   ", (*cont)[i][j]);
                    fprintf(streamc, "\n");
                }
            }
            fclose(streamc);

            FILE *streamt; //ストリームのポインタ設定
            char buffert[30];
            sprintf(buffert, "data/temp/2d%d.csv", istep);
            streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

            for (int i = 0; i <= ndmx; i++)
            {
                fprintf(streamt, "%e   ", (*temp)[i][0]);
                fprintf(streamt, "\n");
            }
            fclose(streamt); //ファイルをクローズ

            printf("%d steps have passed!\n", istep);
            printf("The nominal concnetration is %e\n", c0);
        }

        istep = istep + 1;
        if (istep < nstep)
        {
            goto start;
        }

    end:;
    }
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