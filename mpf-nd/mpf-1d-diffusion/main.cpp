
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
#define ND 100
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;

int nstep = 6401;
int pstep = 400;

double dx = 1.0;
double dtime = 0.02;
double gamma0 = 0.5;
double mobi = 1.0;
double delta = 5.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double S0 = 0.5;

double Dl = 5.0;
double Ds = 0.01;

double temp = 2.0;
double cl = 0.2;

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

double phi[N][ND], phi2[N][ND];
double con[ND], con2[ND], cons[ND], conl[ND];

int phinum;
int phiNum[ND];
int phiIdx[N + 1][ND];

double c0, dc0, cddtt, dev1_s, dev2_s, dev1_l, dev2_l;
double cm0, dcm0, icount, lcount;

int i, j, im, ip, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;

double dF, pddtt, sum1;
double termiikk, termjjkk;

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

    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        if (i <= ND / 2)
        {
            phi[1][i] = 1.0;
            cons[i] = c1e;
            phi[0][i] = 0.0;
            conl[i] = c01e;
        }
        else
        {
            phi[1][i] = 0.0;
            cons[i] = c1e;
            phi[0][i] = 1.0;
            conl[i] = cl;
        }
        con[i] = cons[i] * phi[1][i] + conl[i] * phi[0][i];
        sum1 += con[i];
    }
    c0 = sum1 / ND;
    cm0 = sum1;

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
        cout << "The nominal concnetration is " << c0 << endl;
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
                if (ii == 1 && jj == 0)
                {
                    dF = ((conl[i] - ce) * ml1 + Te - temp) * S0;
                }
                else if (ii == 0 && jj == 1)
                {
                    dF = -(((conl[i] - ce) * ml1 + Te - temp)) * S0;
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

    for (i = 0; i <= ndm; i++)
    {
        if (phi[0][i] == 0.0)
        {
            conl[i] = c01e;
        }
        else if (phi[0][i] > 0.0 && phi[0][i] < 1.0)
        {
            cons[i] = c1e;
            conl[i] = (con[i] - cons[i] * phi[1][i]) / phi[0][i];
            if (conl[i] > 1.0)
            {
                conl[i] = 1.0;
            }
            if (conl[i] < 0.0)
            {
                conl[i] = 0.0;
            }
        }
        else if (phi[0][i] == 1.0)
        {
            conl[i] = con[i];
        }
        con[i] = cons[i] * phi[1][i] + conl[i] * phi[0][i];
        if (con[i] > 1.0)
        {
            con[i] = 1.0;
        }
        if (con[i] < 0.0)
        {
            con[i] = 0.0;
        }
    }

    // Evolution Equation of Concentration field
    for (i = 0; i <= ndm; i++)
    {
        ip = i + 1;
        im = i - 1;
        if (i == ndm)
        {
            ip = ndm;
        }
        if (i == 0)
        {
            im = 0;
        }
        //拡散方程式内における微分計算
        dev1_s = 0.25 * ((phi[1][ip] - phi[1][im]) * (cons[ip] - cons[im])) / dx / dx;
        dev1_l = 0.25 * ((phi[0][ip] - phi[0][im]) * (conl[ip] - conl[im])) / dx / dx;
        dev2_s = phi[1][i] * (cons[ip] + cons[im] - 2.0 * cons[i]) / dx / dx;
        dev2_l = phi[0][i] * (conl[ip] + conl[im] - 2.0 * conl[i]) / dx / dx;

        cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l); //拡散方程式[式(4.42)]
        con2[i] = con[i] + cddtt * dtime;                        //濃度場の時間発展(陽解法)
                                                                 // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
    }

    for (i = 0; i <= ndm; i++)
    {
        con[i] = con2[i]; //補助配列を主配列に移動（濃度場）
    }

    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        sum1 += con[i];
    }
    c0 = sum1 / ND;
    dcm0 = sum1 - cm0;
    icount = 0.0;
    sum1 = 0.0;
    lcount = 0.0;
    // collect mass in interface and liquid
    for (i = 0; i <= ndm; i++)
    {
        if (phi[0][i] > 0.0 && phi[0][i] < 1.0)
        {
            icount += 1.0;
        }
        if (phi[0][i] == 1.0)
        {
            sum1 += con[i];
            lcount += 1.0;
        }
    }
    cl = sum1 / lcount;
    // correction for mass conservation
    for (i = 0; i <= ndm; i++)
    {
        if (phi[0][i] > 0.0 && phi[0][i] < 1.0)
        {
            con[i] = con[i] - dcm0 / icount;
            if (con[i] > 1.0)
            {
                con[i] = 1.0;
            }
            if (con[i] < 0.0)
            {
                con[i] = 0.0;
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
    sprintf(buffer, "data/phi/1d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(stream, "%e   ", phi[1][i]);
        fprintf(stream, "\n");
    }
    fclose(stream); //ファイルをクローズ

    FILE *streamc; //ストリームのポインタ設定
    char bufferc[30];
    sprintf(bufferc, "data/con/1d%d.csv", step);
    streamc = fopen(bufferc, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamc, "%e   ", con[i]);
        fprintf(streamc, "\n");
    }
    fclose(streamc); //ファイルをクローズ

    FILE *streamcl; //ストリームのポインタ設定
    char buffercl[30];
    sprintf(buffercl, "data/conl/1d%d.csv", step);
    streamcl = fopen(buffercl, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamcl, "%e   ", conl[i]);
        fprintf(streamcl, "\n");
    }
    fclose(streamcl); //ファイルをクローズ

    FILE *streamcs; //ストリームのポインタ設定
    char buffercs[30];
    sprintf(buffercs, "data/cons/1d%d.csv", step);
    streamcs = fopen(buffercs, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamcs, "%e   ", cons[i]);
        fprintf(streamcs, "\n");
    }
    fclose(streamcs); //ファイルをクローズ
}
