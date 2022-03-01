
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
#define ND 400
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;

int nstep = 10001;
int pstep = 500;

double dx = 1.0e-7;
double dtime = 4.0e-10;

double gamma0 = 0.1;
double delta = 6.0 * dx;

double mobi = 4.20951e-05;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F0 = 4.0e4;

double As = 1.0e+02;
double cs0 = 0.1;
double Al = 1.0e+01;
double cl0 = 0.4;
double Ds = 0.1e-6;
double Dl = 0.5e-5;

double con[ND][ND], con2[ND][ND], cons[ND][ND], conl[ND][ND];

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

double phi[N][ND][ND], phi2[N][ND][ND];

int phinum;
int phiNum[ND][ND];
int phiIdx[N + 1][ND][ND];

int i, j, im, ip, jm, jp, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;

double pddtt, sum1;
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
    cout << "***************   " << ND * dx << " m"
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
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            if (((i - ND / 2) * (i - ND / 2) + (j - ND / 2) * (j - ND / 2)) <= 16)
            {
                phi[0][i][j] = 1.0;
                cons[i][j] = 0.1;
                phi[1][i][j] = 0.0;
                conl[i][j] = 0.0;
            }
            else
            {
                phi[0][i][j] = 0.0;
                cons[i][j] = 0.0;
                phi[1][i][j] = 1.0;
                conl[i][j] = 0.2;
            }
            con[i][j] = phi[0][i][j] * cons[i][j] + phi[1][i][j] * conl[i][j];
            con2[i][j] = con[i][j];
            sum1 += con[i][j];
        }
    }
    c0 = sum1 / ND / ND;

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            ip = i + 1;
            im = i - 1;
            if (i == ndm)
            {
                ip = 0;
            }
            if (i == 0)
            {
                im = ndm;
            }
            jp = j + 1;
            jm = j - 1;
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
            if (i == ndm)
            {
                ip = 0;
            }
            if (i == 0)
            {
                im = ndm;
            }
            jp = j + 1;
            jm = j - 1;
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

                        termiikk = aij[ii][kk] * (phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j] + phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j]) / (dx * dx);

                        termjjkk = aij[jj][kk] * (phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j] + phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j]) / (dx * dx);

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                    }
                    if (ii != nm && jj == nm)
                    {
                        dF = F0 * (0.4 - conl[i][j]) * 40;
                    }
                    else if (ii == nm && jj != nm)
                    {
                        dF = -F0 * (0.4 - conl[i][j]) * 40;
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

    //
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
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
    for (i = 0; i < ndm; i++)
    {
        for (j = 0; j < ndm; j++)
        {
            cons[i][j] = (Al * con[i][j] + (As * cs0 - Al * cl0) * phi[1][i][j]) / (Al * phi[0][i][j] + As * phi[1][i][j]);
            conl[i][j] = (As * con[i][j] + (Al * cl0 - As * cs0) * phi[0][i][j]) / (Al * phi[0][i][j] + As * phi[1][i][j]);
            if (cons[i][j] >= 1.0)
            {
                cons[i][j] = 1.0;
            }
            if (cons[i][j] <= 0.0)
            {
                cons[i][j] = 0.0;
            }
            if (conl[i][j] >= 1.0)
            {
                conl[i][j] = 1.0;
            }
            if (conl[i][j] <= 0.0)
            {
                conl[i][j] = 0.0;
            }
        }
    }

    // Evolution Equation of Concentration field
    for (i = 1; i <= ndm - 1; i++)
    {
        for (j = 1; j <= ndm - 1; j++)
        {
            ip = i + 1;
            im = i - 1;
            jp = j + 1;
            jm = j - 1;

            //拡散方程式内における微分計算
            dev1_s = 0.25 * ((phi[0][ip][j] - phi[0][im][j]) * (cons[ip][j] - cons[im][j]) + (phi[0][i][jp] - phi[0][i][jm]) * (cons[i][jp] - cons[i][jm])) / dx / dx;
            dev1_l = 0.25 * ((phi[1][ip][j] - phi[1][im][j]) * (conl[ip][j] - conl[im][j]) + (phi[1][i][jp] - phi[1][i][jm]) * (conl[i][jp] - conl[i][jm])) / dx / dx;
            dev2_s = phi[0][i][j] * (cons[ip][j] + cons[im][j] - 2.0 * cons[i][j] + cons[i][jp] + cons[i][jm] - 2.0 * cons[i][j]) / dx / dx;
            dev2_l = phi[1][i][j] * (conl[ip][j] + conl[im][j] - 2.0 * conl[i][j] + conl[i][jp] + conl[i][jm] - 2.0 * conl[i][j]) / dx / dx;

            cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l); //拡散方程式[式(4.42)]
            con2[i][j] = con[i][j] + cddtt * dtime;                  //濃度場の時間発展(陽解法)
                                                                     // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            con[i][j] = con2[i][j]; //補助配列を主配列に移動（濃度場）
        }
    }

    //*** 濃度場の収支補正 *************************************************************
    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            sum1 += con[i][j];
        }
    }
    dc0 = sum1 / ND / ND - c0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            con[i][j] = con[i][j] - dc0;
            if (con[i][j] > 1.0)
            {
                con[i][j] = 1.0;
            }
            if (con[i][j] < 0.0)
            {
                con[i][j] = 0.0;
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

    for (int i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            // double iphi = 0.0;
            // for (int k = 0; k <= nm; k++)
            // {
            //     iphi += phi[k][i] * phi[k][i];
            // }
            fprintf(stream, "%e   ", phi[0][i][j]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream); //ファイルをクローズ

    FILE *streamc;
    char bufferc[30];
    sprintf(bufferc, "data/con/2d%d.csv", step);
    streamc = fopen(bufferc, "a");

    for (int i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            fprintf(streamc, "%e   ", con[i][j]);
            fprintf(streamc, "\n");
        }
    }
    fclose(streamc);
}
