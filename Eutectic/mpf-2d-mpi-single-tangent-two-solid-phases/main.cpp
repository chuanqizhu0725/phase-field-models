
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
#define DRND(x) ((double)(x) / RAND_MAX * rand()) //乱数の関数設定

int nm = N - 1;
int ndm = ND - 1;

int nstep = 1;
int pstep = 1;

double dx = 1.0e-8;
double dtime = 4.0e-12;

double gamma0 = 0.1;
double astre = 0.05;
double delta = 5.0 * dx;

double mobi = 4.2e-5;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F10 = 9.0e6;
double F20 = 2.0e6;

double Ds = 0.1e-6;
double Dl = 0.2e-5;

double temp = 600.0; // K
double Te = 800.0;
double ce = 0.4;
double cl = 0.4;

double ml1 = -2000.0;
double kap1 = 0.25;
double c01e = ce + (temp - Te) / ml1;
double c1e = c01e * kap1;

double ml2 = 2000.0;
double kap2 = 2.5;
double c02e = ce + (temp - Te) / ml2;
double c2e = c02e * kap2;

double c0, dc0;

double phis, phil;

double cddtt, dev1_s, dev2_s, dev1_l, dev2_l;

double con[ND][ND], con_new[ND][ND], con1[ND][ND], con2[ND][ND], con0[ND][ND];

double mij[N][N], aij[N][N], wij[N][N], fij[N][N], anij[N][N], thij[N][N];

double phi[N][ND][ND], phi_new[N][ND][ND];

double phis_ij, phis_ipj, phis_imj, phis_ijp, phis_ijm;
double cons_ij, cons_ipj, cons_imj, cons_ijp, cons_ijm;

int phinum;
int phiNum[ND][ND];
int phiIdx[N + 1][ND][ND];

int i, j, im, ip, jm, jp, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;

double pddtt, sum1;

double theta, theta0;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidxx, phidyy, phidxy;
double ep, ep1p, ep2p;

double dF;

void datasave(int step);
double calcTheta(double dy, double dx);

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
            anij[i][j] = false;
            thij[i][j] = 0.0;
            if ((i == 0) || (j == 0))
            {
                anij[i][j] = true;
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                anij[i][j] = false;
            }
        }
    }

    // Initialize the con field and phase fields
    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            // if (i <= ND / 8 && j < ND / 2)
            if ((i - ND / 2) * (i - ND / 2) + (j - ND / 2) * (j - ND / 2) <= 400)
            {
                phi[1][i][j] = 1.0;
                con1[i][j] = c1e;
                phi[2][i][j] = 0.0;
                con2[i][j] = c2e;
                phi[0][i][j] = 0.0;
                con0[i][j] = c01e;
            }
            // else if (i <= ND / 8 && j >= ND / 2)
            // {
            //     phi[2][i][j] = 1.0;
            //     con2[i][j] = c2e;
            //     phi[1][i][j] = 0.0;
            //     con1[i][j] = c1e;
            //     phi[0][i][j] = 0.0;
            //     con0[i][j] = c02e;
            // }
            else
            {
                phi[1][i][j] = 0.0;
                con1[i][j] = c1e;
                phi[2][i][j] = 0.0;
                con2[i][j] = c2e;
                phi[0][i][j] = 1.0;
                con0[i][j] = cl;
            }
            con[i][j] = phi[1][i][j] * con1[i][j] + phi[2][i][j] * con2[i][j] + phi[0][i][j] * con0[i][j];
            con_new[i][j] = con[i][j];
            sum1 += con[i][j];
        }
    }
    c0 = sum1 / ND / ND;
    cout << "nominal concentration is: " << c0 << endl;

start:;

    // Log in console and output raw data
    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;

        sum1 = 0.0;
        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                sum1 += con[i][j];
            }
        }
        cout << "nominal concentration is: " << sum1 / ND / ND << endl;
    }

    // Calculate local concentration
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            con1[i][j] = c1e;
            con2[i][j] = c2e;
            // liquid in pure phase 1 or phase 2 or their mixer
            if (phi[0][i][j] == 0.0)
            {
                con0[i][j] = c01e * phi[1][i][j] + c02e * phi[2][i][j];
            }
            // liquid in phase 1
            else if (phi[0][i][j] > 0.0 && phi[1][i][j] > 0.0 && phi[2][i][j] == 0.0)
            {
                con0[i][j] = (con[i][j] - phi[1][i][j] * con1[i][j]) / phi[0][i][j];
                if (con0[i][j] >= c01e)
                {
                    con0[i][j] = c01e;
                }
            }
            // liquid in phase 2
            else if (phi[0][i][j] > 0.0 && phi[2][i][j] > 0.0 && phi[1][i][j] == 0.0)
            {
                con0[i][j] = (con[i][j] - phi[2][i][j] * con2[i][j]) / phi[0][i][j];
                if (con0[i][j] <= c02e)
                {
                    con0[i][j] = c02e;
                }
            }
            // liquid in mixer of phase 1 and 2
            else if (phi[0][i][j] > 0.0 && phi[1][i][j] > 0.0 && phi[2][i][j] > 0.0)
            {
                con0[i][j] = (con[i][j] - phi[2][i][j] * con2[i][j] - phi[1][i][j] * con1[i][j]) / phi[0][i][j];
                if (con0[i][j] >= c01e)
                {
                    con0[i][j] = c01e;
                }
                else if (con0[i][j] <= c02e)
                {
                    con0[i][j] = c02e;
                }
            }
            // pure liquid
            else if ((phi[0][i][j] == 1.0 && phi[1][i][j] == 0.0 && phi[2][i][j] == 0.0))
            {
                con0[i][j] = con[i][j];
            }
            // The local concentation should be re-assigned after setting cons and conl (very important)
            con[i][j] = con1[i][j] * phi[1][i][j] + phi[2][i][j] * con2[i][j] + con0[i][j] * phi[0][i][j];
        }
    }

    // Evolution Equation of Concentration field
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
                ip = 0;
            }
            if (i == 0)
            {
                im = ndm;
            }
            if (j == ndm)
            {
                jp = 0;
            }
            if (j == 0)
            {
                jm = ndm;
            }

            phis_ij = phi[1][i][j] + phi[2][i][j];
            phis_ipj = phi[1][ip][j] + phi[2][ip][j];
            phis_imj = phi[1][im][j] + phi[2][im][j];
            phis_ijp = phi[1][i][jp] + phi[2][i][jp];
            phis_ijm = phi[1][i][jm] + phi[2][i][jm];

            cons_ij = con1[i][j] + con2[i][j];
            cons_ipj = con1[ip][j] + con2[ip][j];
            cons_imj = con1[im][j] + con2[im][j];
            cons_ijp = con1[i][jp] + con2[i][jp];
            cons_ijm = con1[i][jm] + con2[i][jm];

            dev1_s = 0.25 * ((phis_ipj - phis_imj) * (cons_ipj - cons_imj) + (phis_ijp - phis_ijm) * (cons_ijp - cons_ijm)) / dx / dx;
            dev1_l = 0.25 * ((phi[0][ip][j] - phi[0][im][j]) * (con0[ip][j] - con0[im][j]) + (phi[0][i][jp] - phi[0][i][jm]) * (con0[i][jp] - con0[i][jm])) / dx / dx;
            dev2_s = (phis_ij) * (cons_ipj + cons_imj + cons_ijp + cons_ijm - 4.0 * cons_ij) / dx / dx;
            dev2_l = phi[0][i][j] * (con0[ip][j] + con0[im][j] + con0[i][jp] + con0[i][jm] - 4.0 * con0[i][j]) / dx / dx;

            cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l);
            con_new[i][j] = con[i][j] + cddtt * dtime;
            // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001;
        }
    }

    // Update con
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            con[i][j] = con_new[i][j];
        }
    }

    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            sum1 += con[i][j];
        }
    }

    // Remove numerical addtion of concentration
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

    // Collection local infomation of phase fields
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
                ip = 0;
            }
            if (i == 0)
            {
                im = ndm;
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
                    ((phi[ii][i][j] == 0.0) && (phi[ii][ip][j] > 0.0) ||
                     (phi[ii][im][j] > 0.0) ||
                     (phi[ii][i][jp] > 0.0) ||
                     (phi[ii][i][jm] > 0.0)))
                {
                    phinum++;
                    phiIdx[phinum][i][j] = ii;
                }
            }
            phiNum[i][j] = phinum;
        }
    }

    FILE *stream;
    char buffer[30];
    sprintf(buffer, "phinum.dat");
    stream = fopen(buffer, "a");

    for (int i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            fprintf(stream, "%d   ", phiNum[i][j]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream);
    FILE *fp;

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
                ip = 0;
            }
            if (i == 0)
            {
                im = ndm;
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

                        phidx = (phi[kk][ip][j] - phi[kk][im][j]) / 2.0 / dx; //フェーズフィールドの空間１階微分
                        phidy = (phi[kk][i][jp] - phi[kk][i][jm]) / 2.0 / dx;
                        phidxx = (phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j]) / (dx * dx); //フェーズフィールドの空間２階微分
                        phidyy = (phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j]) / (dx * dx);
                        phidxy = (phi[kk][ip][jp] + phi[kk][im][jm] - phi[kk][im][jp] - phi[kk][ip][jm]) / 4.0 / (dx * dx);

                        if (anij[ii][kk])
                        {
                            epsilon0 = sqrt(aij[ii][kk]);
                            theta0 = thij[ii][kk];
                            theta = calcTheta(phidy, phidx);                                    //界面の法線方向の角度[式(4.24)]
                            ep = epsilon0 * (1.0 + astre * cos(4.0 * (theta - theta0)));        //勾配エネルギー係数の平方根[式(4.23)]
                            ep1p = -epsilon0 * astre * 4.0 * sin(4.0 * (theta - theta0));       // epの角度による１階微分
                            ep2p = -epsilon0 * astre * 4.0 * 4.0 * cos(4.0 * (theta - theta0)); // epの角度による２階微分

                            termiikk = ep * ep * (phidxx + phidyy) + ep * ep1p * ((phidyy - phidxx) * sin(2.0 * theta) + 2.0 * phidxy * cos(2.0 * theta)) - 0.5 * (ep1p * ep1p + ep * ep2p) * (2.0 * phidxy * sin(2.0 * theta) - phidxx - phidyy - (phidyy - phidxx) * cos(2.0 * theta));
                        }
                        else
                        {
                            termiikk = aij[ii][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]) / (dx * dx);
                        }

                        if (anij[jj][kk])
                        {
                            epsilon0 = sqrt(aij[jj][kk]);
                            theta0 = thij[jj][kk];
                            theta = calcTheta(phidy, phidx);                                    //界面の法線方向の角度[式(4.24)]
                            ep = epsilon0 * (1.0 + astre * cos(4.0 * (theta - theta0)));        //勾配エネルギー係数の平方根[式(4.23)]
                            ep1p = -epsilon0 * astre * 4.0 * sin(4.0 * (theta - theta0));       // epの角度による１階微分
                            ep2p = -epsilon0 * astre * 4.0 * 4.0 * cos(4.0 * (theta - theta0)); // epの角度による２階微分

                            termjjkk = ep * ep * (phidxx + phidyy) + ep * ep1p * ((phidyy - phidxx) * sin(2.0 * theta) + 2.0 * phidxy * cos(2.0 * theta)) - 0.5 * (ep1p * ep1p + ep * ep2p) * (2.0 * phidxy * sin(2.0 * theta) - phidxx - phidyy - (phidyy - phidxx) * cos(2.0 * theta));
                        }
                        else
                        {
                            termjjkk = aij[jj][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]) / (dx * dx);
                        }

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                    }
                    if (ii == 1 && jj == 0)
                    {
                        dF = F10 * (c01e - con0[i][j]) * (1 + 0.2 * (DRND(1.0) - 0.5));
                    }
                    else if (ii == 0 && jj == 1)
                    {
                        dF = -F10 * (c01e - con0[i][j]) * (1 + 0.2 * (DRND(1.0) - 0.5));
                    }
                    else if (ii == 2 && jj == 0)
                    {
                        dF = -F20 * (c02e - con0[i][j]);
                    }
                    else if (ii == 0 && jj == 2)
                    {
                        dF = F20 * (c02e - con0[i][j]);
                    }
                    else
                    {
                        dF = 0.0;
                    }
                    pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j] * phi[jj][i][j]));
                }
                phi_new[ii][i][j] = phi[ii][i][j] + pddtt * dtime;
                if (phi_new[ii][i][j] >= 1.0)
                {
                    phi_new[ii][i][j] = 1.0;
                }
                if (phi_new[ii][i][j] <= 0.0)
                {
                    phi_new[ii][i][j] = 0.0;
                }
            }
        } // i
    }

    // Update phi
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= nm; k++)
            {
                phi[k][i][j] = phi_new[k][i][j];
            }
        }
    }

    // Normalize phi
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
    sprintf(buffer, "data/phi/2d%d.csv", step);
    stream = fopen(buffer, "a");

    for (int i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            fprintf(stream, "%e   ", phi[1][i][j]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream);

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

    FILE *streamcs;
    char buffercs[30];
    sprintf(buffercs, "data/cons/2d%d.csv", step);
    streamcs = fopen(buffercs, "a");

    for (int i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            fprintf(streamcs, "%e   ", con1[i][j]);
            fprintf(streamcs, "\n");
        }
    }
    fclose(streamcs);

    FILE *streamcl;
    char buffercl[30];
    sprintf(buffercl, "data/conl/2d%d.csv", step);
    streamcl = fopen(buffercl, "a");

    for (int i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            fprintf(streamcl, "%e   ", con0[i][j]);
            fprintf(streamcl, "\n");
        }
    }
    fclose(streamcl);
}

double calcTheta(double dy, double dx)
{
    if (dx != 0.0)
    {
        return atan(dy / dx);
    }
    else if (dx == 0.0 && dy > 0)
    {
        return PI / 2.0;
    }
    else if (dx == 0.0 && dy < 0)
    {
        return PI / (-2.0);
    }
    else
    {
        return 0;
    }
}