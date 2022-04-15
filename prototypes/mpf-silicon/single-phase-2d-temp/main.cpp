
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <omp.h>

using namespace std;

#define N 2
#define NDX 80
#define NDY 10
#define NTH 8
#define PI 3.14159

int nm = N - 1;
int ndm = NDX - 1;
int mid = NDX / 4;

int nstep = 500001;
int pstep = 50000;

double dx = 1.0e-5;
double dtime = 1.0e-6;
double gamma0 = 0.5;
double mobi = 1.0e-11;
double delta = 5.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);

double Tm = 1687.0;
double sph_s = 2.29e6;
double kap_s = 22;
double sph_l = 2.53e6;
double kap_l = 54;

double Cdt_s = kap_s / sph_s;
double Cdt_l = kap_l / sph_l;
double dH = 4.122e9;

double Tg = 8.0e3;
double Tv = 0.5e-4;
double Tr = Tg * Tv;

double T_left = 1682.09 - NDX / 4 * dx * Tg;
double T_right = T_left + Tg * NDX * dx;

double mij[N][N], aij[N][N], wij[N][N];

double phi[N][NDX][NDY], phi2[N][NDX][NDY];
double temp[NDX][NDY], temp2[NDX][NDY];
double tempipj, tempimj;

int phinum;
int phiNum[NDX][NDY];
int phiIdx[N + 1][NDX][NDY];

int i, j, im, ip, jp, jm, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;
int intpos, dist, curpos, prepos, frapass;
int curt, pret;
double int_vel;
double int_temp;

double F0, pddtt, sum1;
double termiikk, termjjkk;
double phidxx, phidyy;

double Tddtt;

void datasave(int step);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "temperature field stablity number is: " << dtime * Cdt_s / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * M0 * A0 << endl;

    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
            }
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            if (i <= NDX / 4)
            {
                phi[1][i][j] = 1.0;
                phi[0][i][j] = 0.0;
            }
            else
            {
                phi[1][i][j] = 0.0;
                phi[0][i][j] = 1.0;
            }
            // temp[i][j] = T_left + i * dx * Tg;
        }
    }

#pragma omp parallel num_threads(NTH)
    {

    start:;

        if ((((int)(istep) % pstep) == 0))
        {
            // curpos = frapass + intpos;
            // int_vel = double(curpos - prepos) * dx / double(curt - pret) / dtime;
            // prepos = curpos;
            // pret = curt;
            datasave(istep);
            cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
            // if (int_vel > 0.0)
            // {
            //     cout << "interface velocity is: " << int_vel << endl;
            // }
            // cout << "interface temperature is " << temp[intpos][j] << endl;
        }

        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                if (i == ndm)
                {
                    ip = ndm;
                }
                if (i == 0)
                {
                    im = 0;
                }
                if (j == ndm)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndm;
                }

                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if ((phi[ii][i][j] > 0.0) ||
                        ((phi[ii][i][j] == 0.0) && (phi[ii][ip][j] > 0.0) || (phi[ii][im][j] > 0.0) || (phi[ii][i][jp] > 0.0) || (phi[ii][i][jm] > 0.0)))
                    {
                        phinum++;
                        phiIdx[phinum][i][j] = ii;
                    }
                }
                phiNum[i][j] = phinum;
            }
        }

        // Evolution Equation of Phase Fields
        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                if (i == ndm)
                {
                    ip = ndm;
                }
                if (i == 0)
                {
                    im = 0;
                }
                if (j == ndm)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndm;
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

                            phidxx = (phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j]) / (dx * dx);
                            phidyy = (phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j]) / (dx * dx);

                            termiikk = aij[ii][kk] * phidxx;

                            termjjkk = aij[jj][kk] * phidyy;

                            sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                        }
                        if (ii == 1 && jj == 0)
                        {
                            // F0 = -(temp[i][j] - Tm) * dH / Tm;
                            F0 = 1.0e7;
                        }
                        else if (ii == 0 && jj == 1)
                        {
                            // F0 = (temp[i][j] - Tm) * dH / Tm;
                            F0 = -1.0e7;
                        }
                        pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j]) * (sum1 - 8.0 / PI * F0 * sqrt(phi[ii][i][j] * phi[jj][i][j]));
                    }
                    phi2[ii][i][j] = phi[ii][i][j] + pddtt * dtime;
                    if (phi2[ii][i][j] >= 1.0)
                    {
                        phi2[ii][i][j] = 1.0;
                    }
                    if (phi2[ii][i][j] <= 0.0)
                    {
                        phi2[ii][i][j] = 0.0;
                    }
                    // termperature increase from release of latent heat
                    // if (ii == 1)
                    // {
                    //     temp[i][j] += pddtt * dtime * dH * 1.5 / sph_s;
                    // }
                }
            } // j
        }     // i

        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                sum1 = 0.0;
                for (k = 0; k <= nm; k++)
                {
                    sum1 += phi2[k][i][j];
                }
                for (k = 0; k <= nm; k++)
                {
                    phi2[k][i][j] = phi2[k][i][j] / sum1;
                }
            }
        }

        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                for (k = 0; k <= nm; k++)
                {
                    phi[k][i][j] = phi2[k][i][j];
                }
            }
        }

        // for (i = 0; i <= ndm; i++)
        // {
        //     for (j = 0; j <= ndm; j++)
        //     {
        //         ip = i + 1;
        //         im = i - 1;
        //         jp = j + 1;
        //         jm = j - 1;
        //         // right boundary
        //         if (i == ndm)
        //         {
        //             tempipj = T_right;
        //         }
        //         else
        //         {
        //             tempipj = temp[ip][j];
        //         }
        //         // left boundary
        //         if (i == 0)
        //         {
        //             tempimj = T_left;
        //         }
        //         else
        //         {
        //             tempimj = temp[im][j];
        //         }
        //         if (j == ndm)
        //         {
        //             jp = 0;
        //         }
        //         if (j == 0)
        //         {
        //             jm = ndm;
        //         }
        //         Tddtt = (Cdt_l * phi[0][i][j] + Cdt_s * phi[1][i][j]) * (tempipj + tempimj + temp[i][jp] + temp[i][jm] - 4.0 * temp[i][j]) / dx / dx;
        //         temp2[i][j] = temp[i][j] + Tddtt * dtime;
        //     }
        // }

        // for (i = 0; i <= ndm; i++)
        // {
        //     temp[i][j] = temp2[i][j];
        // }

        // T_left -= Tr * dtime;
        // T_right -= Tr * dtime;

        // intpos = 0;
        // for (i = 0; i <= ndm; i++)
        // {
        //     for (j = 0; j <= ndm; j++)
        //     {
        //         if (phi[1][i][j] < 1.0)
        //         {
        //             intpos = i;
        //             int_temp = temp[i][j];
        //             break;
        //         }
        //     }
        // }

        // if (intpos > mid)
        // {
        //     frapass += 1;
        //     curt = istep;
        //     for (i = 0; i <= ndm - 1; i++)
        //     {
        //         for (j = 0; j <= ndm; j++)
        //         {
        //             for (k = 0; k <= nm; k++)
        //             {
        //                 phi[k][i][j] = phi[k][i + 1][j];
        //             }
        //             temp[i][j] = temp[i + 1][j];
        //         }
        //     }
        //     for (j = 0; j <= ndm; j++)
        //     {
        //         for (k = 0; k <= nm; k++)
        //         {
        //             phi[k][ndm][j] = phi[k][ndm - 1][j];
        //         }
        //         temp[ndm][j] = temp[ndm - 1][j] + Tg * dx;
        //         T_left += Tg * dx;
        //         T_right += Tg * dx;
        //     }
        // }

        istep = istep + 1;
        if (istep < nstep)
        {
            goto start;
        }

    end:;
    }
    return 0;
}

void datasave(int step)
{
    int i, j;
    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/phi/2d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            fprintf(stream, "%e   ", phi[1][i][j]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream); //ファイルをクローズ

    FILE *streamt; //ストリームのポインタ設定
    char buffert[30];
    sprintf(buffert, "data/temp/2d%d.csv", step);
    streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            fprintf(streamt, "%e   ", temp[i][j]);
            fprintf(streamt, "\n");
        }
    }
    fclose(streamt); //ファイルをクローズ
}
