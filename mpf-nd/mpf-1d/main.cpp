
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
#define ND 100
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;

int nstep = 401;
int pstep = 40;

double dx = 1.0;
double dtime = 0.05;
double gamma0 = 0.5;
double mobi = 1.0;
double delta = 7.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F0 = 10.0;

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

void datasave(int step);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    // cout << "concenration field stablity number is: " << dtime * Dl / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * mobi * A0 << endl;
    // cout << "unpper limit of driving force is: " << delta / dtime / mobi / PI << endl;

    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            fij[i][j] = F0;
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                fij[i][j] = 0.0;
            }
            if (i > j)
            {
                fij[i][j] = -fij[i][j];
            }
        }
    }

    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        if (i <= ND / 8)
        {
            phi[0][i] = 1.0;
            phi[1][i] = 0.0;
        }
        else
        {
            phi[0][i] = 0.0;
            phi[1][i] = 1.0;
        }
    }

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
                pddtt += -2.0 * mij[ii][jj] / double(phiNum[i]) * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt(phi[ii][i] * phi[jj][i]));
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
        fprintf(stream, "%e   ", phi[0][i]);
        fprintf(stream, "\n");
    }
    fclose(stream); //ファイルをクローズ
}
