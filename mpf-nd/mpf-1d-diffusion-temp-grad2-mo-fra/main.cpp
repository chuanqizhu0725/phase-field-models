
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

#define N 3
#define ND 100
#define NTH 1
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;
int rows = ND / NTH;

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

double temp0 = 1.0;
double gradT = 0.01;
double cl = 0.2;

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];
double phi[N][ND], phi2[N][ND];
double con[ND], con2[ND], cons[ND], conl[ND];
double temp[ND];
int phiNum[ND];
int phiIdx[N + 1][ND];

void datasave(int step);
double calC01e(double temp0), calC1e(double temp0), calC02e(double temp0), calC2e(double temp0);
double calDF10(double con0, double temp0, double dS), calDF20(double con0, double temp0, double dS);

int main(void)
{

    int i, j, k, im, ip, jm, jp, km, kp;
    double c0, cm0, sum0, lcount, dcm0;

    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "concenration field stablity number is: " << dtime * Dl / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * mobi * A0 << endl;

    // ---------------------------------  Initialization ------------------------------------

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
        temp[i] = temp0 + gradT * i;
    }

    sum0 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        if (i <= ND / 4)
        {
            phi[1][i] = 1.0;
            cons[i] = calC1e(temp[i]);
            phi[0][i] = 0.0;
            conl[i] = calC01e(temp[i]);
        }
        else
        {
            phi[1][i] = 0.0;
            cons[i] = calC1e(temp[i]);
            phi[0][i] = 1.0;
            conl[i] = cl;
        }
        con[i] = cons[i] * phi[1][i] + conl[i] * phi[0][i];
        sum0 += con[i];
    }
    c0 = sum0 / ND;
    cm0 = sum0;

#pragma omp parallel num_threads(NTH)
    {
        int istep, th_id;
        int start, end, offset;
        int ix, ixm, ixp, iy, iym, iyp, ik, ikm, ikp;
        int ii, jj, kk;
        int n1, n2, n3, phinum;

        double cddtt, dev1_s, dev2_s, dev1_l, dev2_l;

        double dF, pddtt, psum, dsum;
        double termiikk, termjjkk;

        th_id = omp_get_thread_num();
        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;
        istep = 0;

    start:;

        // ---------------------------------  Output the calculation data  ------------------------------------

        if (istep % pstep == 0)
        {
            datasave(istep);
            cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
            cout << "The nominal concnetration is " << c0 << endl;
        }

        // ---------------------------------  Collect information of phase fields ------------------------------------
        for (ix = 0; ix <= ndm; ix++)
        {
            ixp = ix + 1;
            ixm = ix - 1;
            if (ix == ndm)
            {
                ixp = ndm;
            }
            if (ix == 0)
            {
                ixm = 0;
            }

            phinum = 0;
            for (ii = 0; ii <= nm; ii++)
            {
                if ((phi[ii][ix] > 0.0) ||
                    ((phi[ii][ix] == 0.0) && (phi[ii][ixp] > 0.0) || (phi[ii][ixm] > 0.0)))
                {
                    phinum++;
                    phiIdx[phinum][ix] = ii;
                }
            }
            phiNum[ix] = phinum;
        }

        // ---------------------------------  Evolution Equation of Phase fields ------------------------------------
        for (ix = 0; ix <= ndm; ix++)
        {
            ixp = ix + 1;
            ixm = ix - 1;
            if (ix == ndm)
            {
                ixp = 0;
            }
            if (ix == 0)
            {
                ixm = ndm;
            }

            for (n1 = 1; n1 <= phiNum[ix]; n1++)
            {
                ii = phiIdx[n1][ix];
                pddtt = 0.0;
                for (n2 = 1; n2 <= phiNum[ix]; n2++)
                {
                    jj = phiIdx[n2][ix];
                    dsum = 0.0;
                    for (n3 = 1; n3 <= phiNum[ix]; n3++)
                    {
                        kk = phiIdx[n3][ix];

                        termiikk = aij[ii][kk] * (phi[kk][ixp] + phi[kk][ixm] - 2.0 * phi[kk][ix]) / (dx * dx);

                        termjjkk = aij[jj][kk] * (phi[kk][ixp] + phi[kk][ixm] - 2.0 * phi[kk][ix]) / (dx * dx);

                        dsum += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][ix];
                    }
                    if (ii == 1 && jj == 0)
                    {
                        dF = calDF10(conl[ix], temp[ix], S0);
                    }
                    else if (ii == 0 && jj == 1)
                    {
                        dF = -calDF10(conl[ix], temp[ix], S0);
                    }
                    else
                    {
                        dF = 0.0;
                    }
                    pddtt += -2.0 * mij[ii][jj] / double(phiNum[ix]) * (dsum - 8.0 / PI * dF * sqrt(phi[ii][ix] * phi[jj][ix]));
                }
                phi2[ii][ix] = phi[ii][ix] + pddtt * dtime;
                if (phi2[ii][ix] >= 1.0)
                {
                    phi2[ii][ix] = 1.0;
                }
                if (phi2[ii][ix] <= 0.0)
                {
                    phi2[ii][ix] = 0.0;
                }
            }
        } // ix

        for (ix = 0; ix <= ndm; ix++)
        {
            for (kk = 0; kk <= nm; kk++)
            {
                phi[kk][ix] = phi2[kk][ix];
            }
        }

        //
        for (ix = 0; ix <= ndm; ix++)
        {
            psum = 0.0;
            for (kk = 0; kk <= nm; kk++)
            {
                psum += phi[kk][ix];
            }
            for (kk = 0; kk <= nm; kk++)
            {
                phi[kk][ix] = phi[kk][ix] / psum;
            }
        }
        // --------------------- Calculate concentration in the interface and liquid phase --------------------------
        for (ix = 0; ix <= ndm; ix++)
        {
            if (phi[0][ix] == 0.0)
            {
                conl[ix] = calC01e(temp[ix]);
            }
            else if (phi[0][ix] > 0.0)
            {
                cons[ix] = calC1e(temp[ix]);
                conl[ix] = (con[ix] - cons[ix] * phi[1][ix]) / phi[0][ix];
                if (conl[ix] > 1.0)
                {
                    conl[ix] = 1.0;
                }
                if (conl[ix] < 0.0)
                {
                    conl[ix] = 0.0;
                }
            }
            con[ix] = cons[ix] * phi[1][ix] + conl[ix] * phi[0][ix];
            if (con[ix] > 1.0)
            {
                con[ix] = 1.0;
            }
            if (con[ix] < 0.0)
            {
                con[ix] = 0.0;
            }
        }

        // --------------------- Correct concentration in liquid phase for mass conservation --------------------------
        sum0 = 0.0;
        lcount = 0.0;
        // collect mass in liquid
        for (ix = 0; ix <= ndm; ix++)
        {
            if (phi[0][ix] == 1.0)
            {
                sum0 += con[ix];
                lcount += phi[0][ix];
            }
        }
        c0 = sum0 / ND;
        dcm0 = sum0 - cm0;
        // correction for mass conservation
        for (ix = 0; ix <= ndm; ix++)
        {
            if (phi[0][ix] == 1.0)
            {
                conl[ix] = conl[ix] - dcm0 / lcount;
                con[ix] = con[ix] - dcm0 / lcount;
                if (con[ix] > 1.0)
                {
                    con[ix] = 1.0;
                }
                if (con[ix] < 0.0)
                {
                    con[ix] = 0.0;
                }
                if (conl[ix] > 1.0)
                {
                    conl[ix] = 1.0;
                }
                if (conl[ix] < 0.0)
                {
                    conl[ix] = 0.0;
                }
            }
        }

        // ---------------------------------  Evolution Equation of Concentration field ------------------------------------
        for (ix = 0; ix <= ndm; ix++)
        {
            ixp = ix + 1;
            ixm = ix - 1;
            if (ix == ndm)
            {
                ixp = ndm;
            }
            if (ix == 0)
            {
                ixm = 0;
            }
            //拡散方程式内における微分計算
            dev1_s = 0.25 * ((phi[1][ixp] - phi[1][ixm]) * (cons[ixp] - cons[ixm])) / dx / dx;
            dev1_l = 0.25 * ((phi[0][ixp] - phi[0][ixm]) * (conl[ixp] - conl[ixm])) / dx / dx;
            dev2_s = phi[1][ix] * (cons[ixp] + cons[ixm] - 2.0 * cons[ix]) / dx / dx;
            dev2_l = phi[0][ix] * (conl[ixp] + conl[ixm] - 2.0 * conl[ix]) / dx / dx;

            cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l); //拡散方程式[式(4.42)]
            con2[ix] = con[ix] + cddtt * dtime;                      //濃度場の時間発展(陽解法)
                                                                     // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
        }

        for (ix = 0; ix <= ndm; ix++)
        {
            con[ix] = con2[ix];
        }

        // moving frame

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

double calC01e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml1 = -10.0;
    double kap1 = 0.2;
    double c01e;
    c01e = ce + (temp0 - Te) / ml1;
    return c01e;
}

double calC1e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml1 = -10.0;
    double kap1 = 0.2;
    double c1e;
    c1e = (ce + (temp0 - Te) / ml1) * kap1;
    return c1e;
}

double calC02e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml2 = 10.0;
    double kap2 = 0.2;
    double c02e;
    c02e = ce + (temp0 - Te) / ml2;
    return c02e;
}

double calC2e(double temp0)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml2 = 10.0;
    double kap2 = 0.2;
    double c2e;
    c2e = 1.0 - (1.0 - (ce + (temp0 - Te) / ml2)) * kap2;
    return c2e;
}

double calDF10(double con0, double temp0, double dS)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml1 = -10.0;
    double DF = ((con0 - ce) * ml1 + Te - temp0) * dS;
    return DF;
}

double calDF20(double con0, double temp0, double dS)
{
    double Te = 0.0;
    double ce = 0.5;
    double ml2 = 10.0;
    double DF = ((con0 - ce) * ml2 + Te - temp0) * dS;
    return DF;
}
