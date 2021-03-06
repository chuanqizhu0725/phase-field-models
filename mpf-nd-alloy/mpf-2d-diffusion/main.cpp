
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

#define N 2
#define NDX 64
#define NDY 64
#define PI 3.14159

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDX - 1;

int nstep = 1001;
int pstep = 200;

double dx = 1.0;
double dtime = 0.02;
double gamma0 = 0.5;
double mobi = 1.0;
double delta = 5.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F0 = 0.5;

double Dl = 5.0;
double Ds = 0.01;

double temp = 2.0;
double cl = 0.2;
// double temp = -1.0;
// double cl = 0.5;

// Linear phae diagram
double Te = 0.0;
double ce = 0.5;
double ml1 = -10.0;
double kap1 = 0.2;
double c01e = ce + (temp - Te) / ml1;
double c1e = c01e * kap1;

double ml2 = 10.0;
double kap2 = 0.2;
double c02e = ce + (temp - Te) / ml2;
double c2e = 1.0 - (1.0 - c02e) * kap2;

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

double phi[N][NDX][NDY], phi2[N][NDX][NDY];
double cont[NDX][NDY], cont2[NDX][NDY], conp[N][NDX][NDY];

int phinum;
int phiNum[NDX][NDY];
int phiIdx[N + 1][NDX][NDY];

double c0, dc0, cddtt, dev01, dev02, dev11, dev12, dev21, dev22;

int i, j, im, ip, jm, jp, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;
int x00, y00, ptype;
double r0, r;

double dF, pddtt, sum1;
double termiikk, termjjkk;
double phidxx, phidyy;

void datasave(int step);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "concenration field stablity number is: " << dtime * Dl / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * mobi * A0 << endl;

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

    // for (i = 0; i <= ndmx; i++)
    // {
    //     for (j = 0; j <= ndmy; j++)
    //     {
    //         phi[0][i][j] = 1.0;
    //         conp[0][i][j] = cl;
    //         for (kk = 1; kk <= nm; kk++)
    //         {

    //             phi[kk][i][j] = 0.0;
    //             conp[kk][i][j] = (kk % 2 == 1) ? c1e : c2e;
    //         }
    //     }
    // }

    // r0 = 10.0;
    // for (ii = 1; ii <= nm; ii++)
    // {
    //     x00 = rand() % NDX;
    //     y00 = rand() % NDY;
    //     for (i = 0; i <= ndmx; i++)
    //     {
    //         for (j = 0; j <= ndmy; j++)
    //         {
    //             r = sqrt(double(((i - x00)) * (i - x00) + (j - y00) * (j - y00)));
    //             if (r <= r0)
    //             {
    //                 if (ii % 2 == 1)
    //                 {
    //                     phi[ii][i][j] = 1.0;
    //                     conp[ii][i][j] = c1e;
    //                     phi[0][i][j] = 0.0;
    //                     conp[0][i][j] = c01e;
    //                     for (kk = 0; kk <= nm; kk++)
    //                     {
    //                         if (kk != ii)
    //                         {
    //                             phi[kk][i][j] = 0.0;
    //                             conp[kk][i][j] = (kk % 2 == 1) ? c1e : c2e;
    //                         }
    //                     }
    //                 }
    //                 else if (ii % 2 == 0)
    //                 {
    //                     phi[ii][i][j] = 1.0;
    //                     conp[ii][i][j] = c2e;
    //                     phi[0][i][j] = 0.0;
    //                     conp[0][i][j] = c02e;
    //                     for (kk = 0; kk <= nm; kk++)
    //                     {
    //                         if (kk != ii)
    //                         {
    //                             phi[kk][i][j] = 0.0;
    //                             conp[kk][i][j] = (kk % 2 == 1) ? c1e : c2e;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // sum1 = 0.0;
    // for (i = 0; i <= ndmx; i++)
    // {
    //     for (j = 0; j <= ndmy; j++)
    //     {
    //         for (kk = 0; kk <= nm; kk++)
    //         {
    //             cont[i][j] += phi[kk][i][j] * conp[kk][i][j];
    //         }
    //         sum1 += cont[i][j];
    //     }
    // }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < 100)
            {
                phi[0][i][j] = 0.0;
                conp[0][i][j] = c01e;
                phi[1][i][j] = 1.0;
                conp[1][i][j] = c1e;
            }
            else
            {
                phi[0][i][j] = 1.0;
                conp[0][i][j] = cl;
                phi[1][i][j] = 0.0;
                conp[1][i][j] = c1e;
            }
            cont[i][j] = conp[1][i][j] * phi[1][i][j] + conp[0][i][j] * phi[0][i][j];
            sum1 += cont[i][j];
        }
    }
    c0 = sum1 / NDX / NDY;

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
        cout << "The nominal concnetration is " << c0 << endl;
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            ip = i + 1;
            im = i - 1;
            jp = j + 1;
            jm = j - 1;
            if (i == ndmx)
            {
                ip = ndmx;
            }
            if (i == 0)
            {
                im = 0;
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
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            ip = i + 1;
            im = i - 1;
            jp = j + 1;
            jm = j - 1;
            if (i == ndmx)
            {
                ip = ndmx;
            }
            if (i == 0)
            {
                im = 0;
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

                        phidxx = (phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j]) / (dx * dx);
                        phidyy = (phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j]) / (dx * dx);

                        termiikk = aij[ii][kk] * phidxx;

                        termjjkk = aij[jj][kk] * phidyy;

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                    }
                    if (ii == 1 && jj == 0)
                    {
                        dF = ((conp[0][i][j] - ce) * ml1 + Te - temp) * 0.5;
                    }
                    else if (ii == 0 && jj == 1)
                    {
                        dF = -((conp[0][i][j] - ce) * ml1 + Te - temp) * 0.5;
                    }
                    else if (ii == 2 && jj == 0)
                    {
                        dF = ((conp[0][i][j] - ce) * ml2 + Te - temp) * 0.5;
                    }
                    else if (ii == 0 && jj == 2)
                    {
                        dF = -((conp[0][i][j] - ce) * ml2 + Te - temp) * 0.5;
                    }
                    else
                    {
                        dF = 0.0;
                    }
                    pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j] * phi[jj][i][j]));
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
            }
        } // i
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= nm; k++)
            {
                phi[k][i][j] = phi2[k][i][j];
            }
        }
    }

    //
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            sum1 = 0.0;
            for (k = 0; k <= nm; k++)
            {
                sum1 += phi[k][i][j];
            }
            for (k = 0; k <= nm; k++)
            {
                phi[k][i][j] = phi[k][i][j] / sum1;
            }
        }
    }

    // Calculate the concentration field in solid and liqud phase
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            conp[1][i][j] = c1e;
            conp[2][i][j] = c2e;
            // liquid in pure phase 1 or phase 2 or their mixer
            if (phi[0][i][j] == 0.0)
            {
                conp[0][i][j] = c01e * phi[1][i][j] + c02e * phi[2][i][j];
            }
            // liquid in phase 2 or phase 1 or their mixer
            else if (phi[0][i][j] > 0.0 || (phi[1][i][j] >= 0.0 || phi[2][i][j] >= 0.0))
            {
                conp[0][i][j] = (cont[i][j] - phi[2][i][j] * conp[2][i][j] - phi[1][i][j] * conp[1][i][j]) / phi[0][i][j];
                if (conp[0][i][j] >= c01e)
                {
                    conp[0][i][j] = c01e;
                }
                else if (conp[0][i][j] <= c02e)
                {
                    conp[0][i][j] = c02e;
                }
            }
            cont[i][j] = conp[2][i][j] * phi[2][i][j] + conp[1][i][j] * phi[1][i][j] + conp[0][i][j] * phi[0][i][j];
        }
    }

    // Evolution Equation of Concentration field
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            ip = i + 1;
            im = i - 1;
            jp = j + 1;
            jm = j - 1;
            if (i == ndmx)
            {
                ip = ndmx;
            }
            if (i == 0)
            {
                im = 0;
            }
            if (j == ndmy)
            {
                jp = 0;
            }
            if (j == 0)
            {
                jm = ndmy;
            }
            dev01 = 0.25 * ((phi[0][ip][j] - phi[0][im][j]) * (conp[0][ip][j] - conp[0][im][j]) + (phi[0][i][jp] - phi[0][i][jm]) * (conp[0][i][jp] - conp[0][i][jm])) / dx / dx;
            dev02 = phi[0][i][j] * (conp[0][ip][j] + conp[0][im][j] + conp[0][i][jp] + conp[0][i][jm] - 4.0 * conp[0][i][j]) / dx / dx;

            dev11 = 0.25 * ((phi[1][ip][j] - phi[1][im][j]) * (conp[1][ip][j] - conp[1][im][j]) + (phi[1][i][jp] - phi[1][i][jm]) * (conp[1][i][jp] - conp[1][i][jm])) / dx / dx;
            dev12 = phi[1][i][j] * (conp[1][ip][j] + conp[1][im][j] + conp[1][i][jp] + conp[1][i][jm] - 4.0 * conp[1][i][j]) / dx / dx;

            dev21 = 0.25 * ((phi[2][ip][j] - phi[2][im][j]) * (conp[2][ip][j] - conp[2][im][j]) + (phi[2][i][jp] - phi[2][i][jm]) * (conp[2][i][jp] - conp[2][i][jm])) / dx / dx;
            dev22 = phi[2][i][j] * (conp[2][ip][j] + conp[2][im][j] + conp[2][i][jp] + conp[2][i][jm] - 4.0 * conp[2][i][j]) / dx / dx;

            cddtt = Ds * (dev11 + dev12) + Ds * (dev21 + dev22) + Dl * (dev01 + dev02);
            cont2[i][j] = cont[i][j] + cddtt * dtime; //????????????????????????(?????????)
                                                      // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //????????????????????????(?????????)
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            cont[i][j] = cont2[i][j]; //????????????????????????????????????????????????
        }
    }

    //*** ???????????????????????? *************************************************************
    sum1 = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            sum1 += cont[i][j];
        }
    }
    c0 = sum1 / NDX / NDY;
    dc0 = sum1 / NDX / NDY - c0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            cont[i][j] = cont[i][j] - dc0;
            if (cont[i][j] > 1.0)
            {
                cont[i][j] = 1.0;
            }
            if (cont[i][j] < 0.0)
            {
                cont[i][j] = 0.0;
            }
        }
    }

    istep = istep + 1;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}

void datasave(int step)
{

    FILE *stream; //????????????????????????????????????
    char buffer[30];
    sprintf(buffer, "data/phi/2d%d.csv", step);
    stream = fopen(buffer, "a"); //????????????????????????????????????????????????????????????

    for (int i = 0; i <= ndmx; i++)
    {
        for (int j = 0; j <= ndmy; j++)
        {
            fprintf(stream, "%e   ", phi[1][i][j]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream); //???????????????????????????

    FILE *streamc;
    char bufferc[30];
    sprintf(bufferc, "data/con/2d%d.csv", step);
    streamc = fopen(bufferc, "a");

    for (int i = 0; i <= ndmx; i++)
    {
        for (int j = 0; j <= ndmy; j++)
        {
            fprintf(streamc, "%e   ", cont[i][j]);
            fprintf(streamc, "\n");
        }
    }
    fclose(streamc);
}
