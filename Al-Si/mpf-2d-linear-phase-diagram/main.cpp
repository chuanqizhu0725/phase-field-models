
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

#define N 3
#define NDX 64
#define NDY 64
#define PI 3.14159

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;

int nstep = 401;
int pstep = 40;

double dx = 1.0e-7;
double dtime = 1.0e-6;

double gamma0 = 0.1;
double delta = 6.0 * dx;

double mobi = 1.20951e-8;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F0 = 5.0e4;

double temp = 880.0; // K
double cl = 0.04;
// double cl = 0.18;

// Linear phae diagram
double Te = 850.0;
double ce = 0.12;
double ml1 = -690.0;
double kap1 = 0.17;
double c01e = ce + (temp - Te) / ml1;
double c1e = c01e * kap1;

double ml2 = 1165.0;
double c02e = ce + (temp - Te) / ml2;
double c2e = 1.0;

double Ds = 1.0e-12;
double Dl = 3.0e-9;

double cont[NDX][NDY], cont2[NDX][NDY], conp[N][NDX][NDY];

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

double phi[N][NDX][NDY], phi2[N][NDX][NDY];

double phis_ij, phis_ipj, phis_imj, phis_ijp, phis_ijm;
double cons_ij, cons_ipj, cons_imj, cons_ijp, cons_ijm;

int phinum;
int phiNum[NDX][NDY];
int phiIdx[N + 1][NDX][NDY];

int i, j, im, ip, jm, jp, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;

double pddtt, sum1;
double phidxx, phidyy;
double termiikk, termjjkk;

double dF;

double c0, dc0, cddtt, dev1_s, dev2_s, dev1_l, dev2_l;

void datasave(int step);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "concenration field stablity number is: " << dtime * Dl / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * mobi * gamma0 << endl;
    cout << "unpper limit of driving force is: " << delta / dtime / mobi / PI << endl;
    cout << "The dimension is" << endl;
    cout << "***************   " << NDX * dx << " m"
         << "    ***************" << endl;
    cout << "**                                          **" << endl;
    cout << "**********************************************" << endl;

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

    sum1 = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            if (i <= NDX / 8)
            {
                phi[1][i][j] = 1.0;
                conp[1][i][j] = c1e;
                phi[2][i][j] = 0.0;
                conp[2][i][j] = c2e;
                phi[0][i][j] = 0.0;
                conp[0][i][j] = c01e;
            }
            else
            {
                phi[1][i][j] = 0.0;
                conp[1][i][j] = c1e;
                phi[2][i][j] = 0.0;
                conp[2][i][j] = c2e;
                phi[0][i][j] = 1.0;
                conp[0][i][j] = cl;
            }
            cont[i][j] = phi[1][i][j] * conp[1][i][j] + phi[2][i][j] * conp[2][i][j] + phi[0][i][j] * conp[0][i][j];
            cont2[i][j] = cont[i][j];
            sum1 += cont[i][j];
        }
    }
    c0 = sum1 / NDX / NDY;
    cout << "nominal concentration is: " << c0 << endl;

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;

        sum1 = 0.0;
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                sum1 += cont[i][j];
            }
        }
        cout << "nominal concentration is: " << sum1 / NDX / NDY << endl;
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
                        dF = F0 * ((conp[0][i][j] - ce) * ml1 + Te - temp);
                    }
                    else if (ii == 0 && jj == 1)
                    {
                        dF = -F0 * ((conp[0][i][j] - ce) * ml1 + Te - temp);
                    }
                    else if (ii == 2 && jj == 0)
                    {
                        dF = F0 * ((conp[0][i][j] - ce) * ml2 + Te - temp);
                    }
                    else if (ii == 0 && jj == 2)
                    {
                        dF = -F0 * ((conp[0][i][j] - ce) * ml2 + Te - temp);
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
            // liquid
            if (phi[0][i][j] == 1.0)
            {
                conp[1][i][j] = c1e;
                conp[2][i][j] = c2e;
                conp[0][i][j] = cont[i][j];
            }
            // interface of phase 1
            else if (phi[0][i][j] > 0.0 && phi[0][i][j] < 1.0 && phi[1][i][j] > 0.0 && phi[1][i][j] < 1.0 && phi[2][i][j] == 0.0)
            {
                conp[1][i][j] = c1e;
                conp[0][i][j] = (cont[i][j] - conp[1][i][j] * phi[1][i][j]) / phi[0][i][j];
                if (conp[0][i][j] > c01e)
                {
                    conp[0][i][j] = c01e;
                }
            }
            // interface of phase 2
            else if (phi[0][i][j] > 0.0 && phi[0][i][j] < 1.0 && phi[2][i][j] > 0.0 && phi[2][i][j] < 1.0 && phi[1][i][j] == 0.0)
            {
                conp[2][i][j] = c2e;
                conp[0][i][j] = (cont[i][j] - conp[2][i][j] * phi[2][i][j]) / phi[0][i][j];
                if (conp[0][i][j] < c02e)
                {
                    conp[0][i][j] = c02e;
                }
            }
            // phase 1
            else if (phi[1][i][j] == 1.0)
            {
                conp[1][i][j] = c1e;
                conp[0][i][j] = c01e;
            }
            // phase 2
            else if (phi[2][i][j] == 1.0)
            {
                conp[2][i][j] = c2e;
                conp[0][i][j] = c02e;
            }
            cont[i][j] = phi[0][i][j] * conp[0][i][j] + phi[1][i][j] * conp[1][i][j] + phi[2][i][j] * conp[2][i][j];
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
            //拡散方程式内における微分計算

            phis_ij = phi[1][i][j] + phi[2][i][j];
            phis_ipj = phi[1][ip][j] + phi[2][ip][j];
            phis_imj = phi[1][im][j] + phi[2][im][j];
            phis_ijp = phi[1][i][jp] + phi[2][i][jp];
            phis_ijm = phi[1][i][jm] + phi[2][i][jm];

            cons_ij = conp[1][i][j] + conp[2][i][j];
            cons_ipj = conp[1][ip][j] + conp[2][ip][j];
            cons_imj = conp[1][im][j] + conp[2][im][j];
            cons_ijp = conp[1][i][jp] + conp[2][i][jp];
            cons_ijm = conp[1][i][jm] + conp[2][i][jm];

            dev1_s = 0.25 * ((phis_ipj - phis_imj) * (cons_ipj - cons_imj) + (phis_ijp - phis_ijm) * (cons_ijp - cons_ijm)) / dx / dx;
            dev1_l = 0.25 * ((phi[0][ip][j] - phi[0][im][j]) * (conp[0][ip][j] - conp[0][im][j]) + (phi[0][i][jp] - phi[0][i][jm]) * (conp[0][i][jp] - conp[0][i][jm])) / dx / dx;
            dev2_s = phis_ij * (cons_ipj + cons_imj + cons_ijp + cons_ijm - 4.0 * cons_ij) / dx / dx;
            dev2_l = phi[0][i][j] * (conp[0][ip][j] + conp[0][im][j] + conp[0][i][jp] + conp[0][i][jm] - 4.0 * conp[0][i][j]) / dx / dx;

            cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l); //拡散方程式[式(4.42)]
            cont2[i][j] = cont[i][j] + cddtt * dtime;                //濃度場の時間発展(陽解法)
                                                                     // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            cont[i][j] = cont2[i][j]; //補助配列を主配列に移動（濃度場）
        }
    }

    //*** 濃度場の収支補正 *************************************************************
    sum1 = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            sum1 += cont[i][j];
        }
    }
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
    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/phi/2d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndmx; i++)
    {
        for (int j = 0; j <= ndmy; j++)
        {
            fprintf(stream, "%e   ", phi[1][i][j]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream); //ファイルをクローズ

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
