
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
#define ND 400
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;

int nstep = 50001;
int pstep = 10000;

double dx = 1.0e-7;
double dtime = 4.0e-10;

double gamma0 = 0.1;
double delta = 5.0 * dx;

double mobi = 4.20951e-05;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F0 = 4.0e5;

double Ds = 0.1e-6;
double Dl = 0.5e-5;

double temp = 600.0; // K
double Te = 800.0;
double ce = 0.4;

double ml1 = -2000.0;
double kap1 = 0.25;
double c010 = ce + (temp - Te) / ml1;
double c10 = c010 * kap1;

double ml2 = 2000.0;
double kap2 = 2.5;
double c020 = ce + (temp - Te) / ml2;
double c20 = c020 * kap2;

double con[ND], con_new[ND], con1[ND], con2[ND], con0[ND];

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

double phi[N][ND], phi2[N][ND];

int phinum;
int phiNum[ND];
int phiIdx[N + 1][ND];

int i, j, im, ip, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;
double dd;

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
            con1[i] = c10;
            phi[2][i] = 0.0;
            con2[i] = c20;
            phi[0][i] = 0.0;
            con0[i] = c010;
        }
        // else if (i >= ND * 7 / 8)
        // {
        //     phi[2][i] = 1.0;
        //     con2[i] = c20;
        //     phi[1][i] = 0.0;
        //     con1[i] = c10;
        //     phi[0][i] = 0.0;
        //     con0[i] = c020;
        // }
        else
        {
            phi[1][i] = 0.0;
            con1[i] = c10;
            phi[2][i] = 0.0;
            con2[i] = c20;
            phi[0][i] = 1.0;
            con0[i] = 0.2;
        }
        con[i] = phi[1][i] * con1[i] + phi[2][i] * con2[i] + phi[0][i] * con0[i];
        con_new[i] = con[i];
        sum1 += con[i];
    }
    c0 = sum1 / ND;
    cout << "nominal concentration is: " << c0 << endl;

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

    // Calculate the concentration field in solid and liqud phase based on local concentration and equilibrium rule.
    for (i = 0; i <= ndm; i++)
    {
        // pure solid and solid in mixer
        con1[i] = c10;
        con2[i] = c20;
        // liquid in pure phase 1 or phase 2 or their mixer
        if (phi[0][i] == 0.0)
        {
            con0[i] = c010 * phi[1][0] + c020 * phi[2][0];
        }
        // liquid in a mixer with phase 1
        else if (phi[0][i] > 0.0 && phi[2][i] == 0.0)
        {
            con0[i] = (con[i] - phi[1][i] * con1[i]) / phi[0][i];
            if (con0[i] >= c010)
            {
                con0[i] = c010;
            }
        }
        // liquid in a mixer with phase 2
        else if (phi[0][i] > 0.0 && phi[1][i] == 0.0)
        {
            con0[i] = (con[i] - phi[2][i] * con2[i]) / phi[0][i];
            if (con0[i] <= c020)
            {
                con0[i] = c020;
            }
        }
        // liquid in a mixer with phase 2 and phase 1
        else if (phi[0][i] > 0.0 && phi[1][i] >= 0.0 && phi[2][i] >= 0.0)
        {
            con0[i] = (con[i] - phi[2][i] * con2[i] - phi[1][i] * con1[i]) / phi[0][i];
            if (con0[i] > 1.0)
            {
                con0[i] = c0;
            }
        }
        // The local concentation should be re-assigned after setting cons and conl (very important)
        con[i] = con1[i] * phi[1][i] + phi[2][i] * con2[i] + con0[i] * phi[0][i];
    }

    // Evolution Equation of Concentration field
    for (i = 1; i <= ndm - 1; i++)
    {
        ip = i + 1;
        im = i - 1;
        dev1_s = 0.25 * ((phi[1][ip] + phi[2][ip] - phi[1][im] - phi[2][im]) * (con1[ip] + con2[ip] - con1[im] - con2[im])) / dx / dx;
        dev1_l = 0.25 * ((phi[0][ip] - phi[0][im]) * (con0[ip] - con0[im])) / dx / dx;
        dev2_s = (phi[1][i] + phi[2][i]) * (con1[ip] + con2[ip] + con1[im] + con2[im] - 2.0 * con1[i] - 2.0 * con2[i]) / dx / dx;
        dev2_l = phi[0][i] * (con0[ip] + con0[im] - 2.0 * con0[i]) / dx / dx;

        cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l);
        con_new[i] = con[i] + cddtt * dtime;
        // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001;
    }

    // Set the boundary value to a value of inner grid (isolated box)
    con_new[ndm] = con_new[ndm - 1];
    con_new[0] = con_new[1];
    for (i = 0; i <= ndm; i++)
    {
        con[i] = con_new[i];
    }

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
                    dF = F0 * (c010 - con0[i]);
                }
                else if (ii == 0 && jj == 1)
                {
                    dF = -F0 * (c010 - con0[i]);
                }
                else if (ii == 2 && jj == 0)
                {
                    dF = -F0 * (c020 - con0[i]);
                }
                else if (ii == 0 && jj == 2)
                {
                    dF = F0 * (c020 - con0[i]);
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
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "data/phi/1d%d.csv", step);
    stream = fopen(buffer, "a");

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(stream, "%e   ", phi[1][i]);
        fprintf(stream, "\n");
    }
    fclose(stream);

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
        fprintf(streamcs, "%e   ", con1[i]);
        fprintf(streamcs, "\n");
    }
    fclose(streamcs);

    FILE *streamcl;
    char buffercl[30];
    sprintf(buffercl, "data/conl/1d%d.csv", step);
    streamcl = fopen(buffercl, "a");

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamcl, "%e   ", con0[i]);
        fprintf(streamcl, "\n");
    }
    fclose(streamcl);
}
