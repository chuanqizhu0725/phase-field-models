
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
#define ND 256
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;

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

double cont[ND], cont2[ND], conp[N][ND];

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

double phi[N][ND], phi2[N][ND];

double phis_i, phis_ip, phis_im;
double cons_i, cons_ip, cons_im;

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
            phi[1][i] = 1.0;
            conp[1][i] = c1e;
            phi[2][i] = 0.0;
            conp[2][i] = c2e;
            phi[0][i] = 0.0;
            conp[0][i] = c01e;
        }
        else
        {
            phi[1][i] = 0.0;
            conp[1][i] = c1e;
            phi[2][i] = 0.0;
            conp[2][i] = c2e;
            phi[0][i] = 1.0;
            conp[0][i] = cl;
        }
        cont[i] = phi[1][i] * conp[1][i] + phi[2][i] * conp[2][i] + phi[0][i] * conp[0][i];
        cont2[i] = cont[i];
        sum1 += cont[i];
    }
    c0 = sum1 / ND;
    cout << "nominal concentration is: " << c0 << endl;

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;

        sum1 = 0.0;
        for (i = 0; i <= ndm; i++)
        {
            sum1 += cont[i];
        }
        cout << "nominal concentration is: " << sum1 / ND << endl;
    }

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
            ip = ndm;
        }
        if (i == 0)
        {
            im = 0;
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
                    dF = F0 * ((conp[0][i] - ce) * ml1 + Te - temp);
                }
                else if (ii == 0 && jj == 1)
                {
                    dF = -F0 * ((conp[0][i] - ce) * ml1 + Te - temp);
                }
                else if (ii == 2 && jj == 0)
                {
                    dF = F0 * ((conp[0][i] - ce) * ml2 + Te - temp);
                }
                else if (ii == 0 && jj == 2)
                {
                    dF = -F0 * ((conp[0][i] - ce) * ml2 + Te - temp);
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
    for (i = 0; i <= ndm; i++)
    {
        // liquid
        if (phi[0][i] == 1.0)
        {
            conp[1][i] = c1e;
            conp[2][i] = c2e;
            conp[0][i] = cont[i];
        }
        // interface of phase 1
        else if (phi[0][i] > 0.0 && phi[0][i] < 1.0 && phi[1][i] > 0.0 && phi[1][i] < 1.0 && phi[2][i] == 0.0)
        {
            conp[1][i] = c1e;
            conp[0][i] = (cont[i] - conp[1][i] * phi[1][i]) / phi[0][i];
            if (conp[0][i] > c01e)
            {
                conp[0][i] = c01e;
            }
        }
        // interface of phase 2
        else if (phi[0][i] > 0.0 && phi[0][i] < 1.0 && phi[2][i] > 0.0 && phi[2][i] < 1.0 && phi[1][i] == 0.0)
        {
            conp[2][i] = c2e;
            conp[0][i] = (cont[i] - conp[2][i] * phi[2][i]) / phi[0][i];
            if (conp[0][i] < c02e)
            {
                conp[0][i] = c02e;
            }
        }
        // phase 1
        else if (phi[1][i] == 1.0)
        {
            conp[1][i] = c1e;
            conp[0][i] = c01e;
        }
        // phase 2
        else if (phi[2][i] == 1.0)
        {
            conp[2][i] = c2e;
            conp[0][i] = c02e;
        }
        cont[i] = phi[0][i] * conp[0][i] + phi[1][i] * conp[1][i] + phi[2][i] * conp[2][i];
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

        phis_i = phi[1][i] + phi[2][i];
        phis_ip = phi[1][ip] + phi[2][ip];
        phis_im = phi[1][im] + phi[2][im];

        cons_i = conp[1][i] + conp[2][i];
        cons_ip = conp[1][ip] + conp[2][ip];
        cons_im = conp[1][im] + conp[2][im];

        dev1_s = 0.25 * ((phis_ip - phis_im) * (cons_ip - cons_im)) / dx / dx;
        dev1_l = 0.25 * ((phi[0][ip] - phi[0][im]) * (conp[0][ip] - conp[0][im])) / dx / dx;
        dev2_s = phis_i * (cons_ip + cons_im - 2.0 * cons_i) / dx / dx;
        dev2_l = phi[0][i] * (conp[0][ip] + conp[0][im] - 2.0 * conp[0][i]) / dx / dx;

        cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l); //拡散方程式[式(4.42)]
        cont2[i] = cont[i] + cddtt * dtime;                      //濃度場の時間発展(陽解法)
                                                                 // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
    }

    for (i = 0; i <= ndm; i++)
    {
        cont[i] = cont2[i]; //補助配列を主配列に移動（濃度場）
    }

    //*** 濃度場の収支補正 *************************************************************
    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        sum1 += cont[i];
    }
    dc0 = sum1 / ND - c0;
    for (i = 0; i <= ndm; i++)
    {
        cont[i] = cont[i] - dc0;
        if (cont[i] > 1.0)
        {
            cont[i] = 1.0;
        }
        if (cont[i] < 0.0)
        {
            cont[i] = 0.0;
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

    FILE *streamc;
    char bufferc[30];
    sprintf(bufferc, "data/con/1d%d.csv", step);
    streamc = fopen(bufferc, "a");

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamc, "%e   ", cont[i]);
        fprintf(streamc, "\n");
    }
    fclose(streamc);
}
