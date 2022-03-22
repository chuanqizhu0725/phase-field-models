
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
#define ND 800
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;

int nstep = 400001;
int pstep = 10000;

double dx = 1.0e-7;
double dtime = 4.0e-10;

double gamma0 = 0.1;
double delta = 6.0 * dx;

double mobi = 4.20951e-05;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F0 = 5.0e4;

double As = 1.0e+02;
double cs0 = 0.1;
double Al = 1.0e+01;
double cl0 = 0.4;
double Ds = 0.1e-6;
double Dl = 0.5e-5;

double con[ND], con2[ND], cons[ND], conl[ND];

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

double phi[N][ND], phi2[N][ND];

int phinum;
int phiNum[ND];
int phiIdx[N + 1][ND];

int i, j, im, ip, k;
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
        if (i <= ND / 8)
        {
            phi[0][i] = 1.0;
            cons[i] = 0.1;
            phi[1][i] = 0.0;
            conl[i] = 0.0;
        }
        else
        {
            phi[0][i] = 0.0;
            cons[i] = 0.0;
            phi[1][i] = 1.0;
            conl[i] = 0.2;
        }
        con[i] = phi[0][i] * cons[i] + phi[1][i] * conl[i];
        con2[i] = con[i];
        sum1 += con[i];
    }
    c0 = sum1 / ND;

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
    }

    for (i = 0; i <= ndm; i++)
    {
        ip = i + 1;
        im = i - 1;
        if (i == ndm)
        {
            ip = ndm - 1;
        }
        if (i == 0)
        {
            im = 1;
        }

        phinum = 0;
        for (ii = 0; ii <= nm; ii++)
        {
            if ((phi[ii][i] > 0.0) ||
                ((phi[ii][i] == 0.0) && (phi[ii][ip] > 0.0) || (phi[ii][im] > 0.0)))
            {
                phinum++;
                phiIdx[phinum][i] = ii;
            }
        }
        phiNum[i] = phinum;
    }

    // Evolution Equation of Phase Fields
    for (i = 0; i <= ndm; i++)
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

        for (n1 = 1; n1 <= phiNum[i]; n1++)
        {
            ii = phiIdx[n1][i];
            pddtt = 0.0;
            for (n2 = 1; n2 <= phiNum[i]; n2++)
            {
                jj = phiIdx[n2][i];
                sum1 = 0.0;
                for (n3 = 1; n3 <= phiNum[i]; n3++)
                {
                    kk = phiIdx[n3][i];

                    termiikk = aij[ii][kk] * (phi[kk][ip] + phi[kk][im] - 2.0 * phi[kk][i]) / (dx * dx);

                    termjjkk = aij[jj][kk] * (phi[kk][ip] + phi[kk][im] - 2.0 * phi[kk][i]) / (dx * dx);

                    sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i];
                }
                if (ii != nm && jj == nm)
                {
                    dF = F0 * (0.4 - conl[i]) * 40;
                }
                else if (ii == nm && jj != nm)
                {
                    dF = -F0 * (0.4 - conl[i]) * 40;
                }
                else
                {
                    dF = 0.0;
                }
                pddtt += -2.0 * mij[ii][jj] / double(phiNum[i]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i] * phi[jj][i]));
            }
            phi2[ii][i] = phi[ii][i] + pddtt * dtime;
            if (phi2[ii][i] >= 1.0)
            {
                phi2[ii][i] = 1.0;
            }
            if (phi2[ii][i] <= 0.0)
            {
                phi2[ii][i] = 0.0;
            }
        }
    } // i

    for (i = 0; i <= ndm; i++)
    {
        for (k = 0; k <= nm; k++)
        {
            phi[k][i] = phi2[k][i];
        }
    }

    //
    for (i = 0; i <= ndm; i++)
    {
        sum1 = 0.0;
        for (k = 0; k <= nm; k++)
        {
            sum1 += phi[k][i];
        }
        for (k = 0; k <= nm; k++)
        {
            phi[k][i] = phi[k][i] / sum1;
        }
    }

    // Calculate the concentration field in solid and liqud phase
    for (i = 0; i < ndm; i++)
    {
        cons[i] = (Al * con[i] + (As * cs0 - Al * cl0) * phi[1][i]) / (Al * phi[0][i] + As * phi[1][i]);
        conl[i] = (As * con[i] + (Al * cl0 - As * cs0) * phi[0][i]) / (Al * phi[0][i] + As * phi[1][i]);
        if (cons[i] >= 1.0)
        {
            cons[i] = 1.0;
        }
        if (cons[i] <= 0.0)
        {
            cons[i] = 0.0;
        }
        if (conl[i] >= 1.0)
        {
            conl[i] = 1.0;
        }
        if (conl[i] <= 0.0)
        {
            conl[i] = 0.0;
        }
    }

    // Evolution Equation of Concentration field
    for (i = 1; i <= ndm - 1; i++)
    {
        ip = i + 1;
        im = i - 1;
        //拡散方程式内における微分計算
        dev1_s = 0.25 * ((phi[0][ip] - phi[0][im]) * (cons[ip] - cons[im])) / dx / dx;
        dev1_l = 0.25 * ((phi[1][ip] - phi[1][im]) * (conl[ip] - conl[im])) / dx / dx;
        dev2_s = phi[0][i] * (cons[ip] + cons[im] - 2.0 * cons[i]) / dx / dx;
        dev2_l = phi[1][i] * (conl[ip] + conl[im] - 2.0 * conl[i]) / dx / dx;

        cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l); //拡散方程式[式(4.42)]
        con2[i] = con[i] + cddtt * dtime;                        //濃度場の時間発展(陽解法)
                                                                 // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
    }

    for (i = 0; i <= ndm; i++)
    {
        con[i] = con2[i]; //補助配列を主配列に移動（濃度場）
    }

    //*** 濃度場の収支補正 *************************************************************
    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        sum1 += con[i];
    }
    dc0 = sum1 / ND - c0;
    for (i = 0; i <= ndm; i++)
    {
        con[i] = con[i] - dc0;
        if (con[i] > 1.0)
        {
            con[i] = 1.0;
        }
        if (con[i] < 0.0)
        {
            con[i] = 0.0;
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
    sprintf(buffer, "data/phi/1d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        // double iphi = 0.0;
        // for (int k = 0; k <= nm; k++)
        // {
        //     iphi += phi[k][i] * phi[k][i];
        // }
        fprintf(stream, "%e   ", phi[0][i]);
        fprintf(stream, "\n");
    }
    fclose(stream); //ファイルをクローズ

    FILE *streamc;
    char bufferc[30];
    sprintf(bufferc, "data/con/1d%d.csv", step);
    streamc = fopen(bufferc, "a");

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamc, "%e   ", con[i]);
        fprintf(streamc, "\n");
    }
    fclose(streamc);

    FILE *streamcs;
    char buffercs[30];
    sprintf(buffercs, "data/cons/1d%d.csv", step);
    streamcs = fopen(buffercs, "a");

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamcs, "%e   ", cons[i]);
        fprintf(streamcs, "\n");
    }
    fclose(streamcs);

    FILE *streamcl;
    char buffercl[30];
    sprintf(buffercl, "data/conl/1d%d.csv", step);
    streamcl = fopen(buffercl, "a");

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamcl, "%e   ", conl[i]);
        fprintf(streamcl, "\n");
    }
    fclose(streamcl);
}
