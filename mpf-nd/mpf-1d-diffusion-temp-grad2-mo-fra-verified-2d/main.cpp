
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <omp.h>

using namespace std;

#define N 3
#define NTH 4
#define NDX 128
#define NDY 128
#define NDL 2560
#define PI 3.14159

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndml = NDL - 1;
int rows = NDX / NTH;
int rowsl = NDL / NTH;

int nstep = 10001;
int pstep = 1000;

double dx = 1.0;
double dtime = 0.02;
double gamma0 = 0.5;
double mobi = 1.0;
double delta = 5.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double S0 = 0.5;

double Dl = 1.0;
double Ds = 0.001;

double temp0 = 2.0;
double gradT = 0.00;
double dts = 0.00 / nstep;
double cl = 0.8;

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];
double phi[N][NDX][NDY], phi2[N][NDX][NDY];
double cont[NDX][NDY], cont2[NDX][NDY], conp[N][NDX][NDY];
double conlr[NDL], conlr2[NDL];
double temp[NDX][NDY];
int phiNum[NDX][NDY];
int phiIdx[N + 1][NDX][NDY];

void datasave(int step);
double calC01e(double temp0), calC1e(double temp0), calC02e(double temp0), calC2e(double temp0);
double calDF10(double con0, double temp0, double dS), calDF20(double con0, double temp0, double dS);

int main(void)
{

    int i, j, k, im, ip, jm, jp, km, kp;
    int ni, phinum0;
    int intpos, dist, hasS, allS, allL;
    double c0, c00, dc0, sum0;

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

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            temp[i][j] = temp0 + gradT * i;
        }
    }

    sum0 = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            if (i <= NDX / 16)
            {
                phi[1][i][j] = 0.0;
                conp[1][i][j] = calC1e(temp[i][j]);
                phi[2][i][j] = 1.0;
                conp[2][i][j] = calC2e(temp[i][j]);
                phi[0][i][j] = 0.0;
                conp[0][i][j] = calC02e(temp[i][j]);
            }
            else
            {
                phi[1][i][j] = 0.0;
                conp[1][i][j] = calC1e(temp[i][j]);
                phi[2][i][j] = 0.0;
                conp[2][i][j] = calC2e(temp[i][j]);
                phi[0][i][j] = 1.0;
                conp[0][i][j] = cl;
            }
            cont[i][j] = conp[1][i][j] * phi[1][i][j] + conp[2][i][j] * phi[2][i][j] + conp[0][i][j] * phi[0][i][j];
            sum0 += cont[i][j];
        }
    }
    c0 = sum0 / NDX / NDY;
    // Long range field
    // for (i = 0; i <= ndml; i++)
    // {
    //     conlr[i] = cl;
    // }
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

            phinum0 = 0;
            for (ni = 0; ni <= nm; ni++)
            {
                if ((phi[ni][i][j] > 0.0) ||
                    ((phi[ni][i][j] == 0.0) && (phi[ni][ip][j] > 0.0) || (phi[ni][im][j] > 0.0) || (phi[ni][i][jp] > 0.0) || (phi[ni][i][jm] > 0.0)))
                {
                    phinum0++;
                    phiIdx[phinum0][i][j] = ni;
                }
            }
            phiNum[i][j] = phinum0;
        }
    }

#pragma omp parallel num_threads(NTH)
    {
        int istep, th_id;
        int start, end, offset;
        int startl, endL, offsetl;
        int ix, ixm, ixp, iy, iym, iyp, ik, ikm, ikp;
        int ii, jj, kk;
        int n1, n2, n3, phinum;

        double cddtt, sumcs, sumcl;

        double dF, pddtt, psum, dsum;
        double termiikk, termjjkk;

        istep = 0;
        th_id = omp_get_thread_num();

        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;

        offsetl = th_id * rowsl;
        startl = offsetl;
        endL = offsetl + rowsl - 1;

    start:;

        // ---------------------------------  Output the calculation data  ------------------------------------

        if ((istep % pstep == 0) && (th_id == 0))
        {
            datasave(istep);
            cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
            cout << "The nominal concnetration is " << c00 << endl;
        }

        // ---------------------------------  Evolution Equation of Phase fields ------------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                if (ix == ndmx)
                {
                    ixp = 0;
                }
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (iy == ndmy)
                {
                    iyp = 0;
                }
                if (iy == 0)
                {
                    iym = ndmy;
                }

                for (n1 = 1; n1 <= phiNum[ix][iy]; n1++)
                {
                    ii = phiIdx[n1][ix][iy];
                    pddtt = 0.0;
                    for (n2 = 1; n2 <= phiNum[ix][iy]; n2++)
                    {
                        jj = phiIdx[n2][ix][iy];
                        dsum = 0.0;
                        for (n3 = 1; n3 <= phiNum[ix][iy]; n3++)
                        {
                            kk = phiIdx[n3][ix][iy];

                            termiikk = aij[ii][kk] * (phi[kk][ixp][iy] + phi[kk][ixm][iy] + phi[kk][ix][iyp] + phi[kk][ix][iym] - 4.0 * phi[kk][ix][iy]) / (dx * dx);

                            termjjkk = aij[jj][kk] * (phi[kk][ixp][iy] + phi[kk][ixm][iy] + phi[kk][ix][iyp] + phi[kk][ix][iym] - 4.0 * phi[kk][ix][iy]) / (dx * dx);

                            dsum += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][ix][iy];
                        }
                        if (ii == 1 && jj == 0)
                        {
                            dF = calDF10(conp[0][ix][iy], temp[ix][iy], S0);
                        }
                        else if (ii == 0 && jj == 1)
                        {
                            dF = -calDF10(conp[0][ix][iy], temp[ix][iy], S0);
                        }
                        else if (ii == 2 && jj == 0)
                        {
                            dF = calDF20(conp[0][ix][iy], temp[ix][iy], S0);
                        }
                        else if (ii == 0 && jj == 2)
                        {
                            dF = -calDF20(conp[0][ix][iy], temp[ix][iy], S0);
                        }
                        else
                        {
                            dF = 0.0;
                        }
                        pddtt += -2.0 * mij[ii][jj] / double(phiNum[ix][iy]) * (dsum - 8.0 / PI * dF * sqrt(phi[ii][ix][iy] * phi[jj][ix][iy]));
                    }
                    phi2[ii][ix][iy] = phi[ii][ix][iy] + pddtt * dtime;
                    if (phi2[ii][ix][iy] >= 1.0)
                    {
                        phi2[ii][ix][iy] = 1.0;
                    }
                    if (phi2[ii][ix][iy] <= 0.0)
                    {
                        phi2[ii][ix][iy] = 0.0;
                    }
                }
            } // ix
        }     // iy

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (kk = 0; kk <= nm; kk++)
                {
                    phi[kk][ix][iy] = phi2[kk][ix][iy];
                }
            }
        }

        //
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                psum = 0.0;
                for (kk = 0; kk <= nm; kk++)
                {
                    psum += phi[kk][ix][iy];
                }
                for (kk = 0; kk <= nm; kk++)
                {
                    phi[kk][ix][iy] = phi[kk][ix][iy] / psum;
                    phi2[kk][ix][iy] = phi[kk][ix][iy];
                }
            }
        }
#pragma omp barrier
        // ---------------------------  Collect information of phase fields ----------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                if (ix == ndmx)
                {
                    ixp = 0;
                }
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (iy == ndmy)
                {
                    iyp = 0;
                }
                if (iy == 0)
                {
                    iym = ndmy;
                }

                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if ((phi[ii][ix][iy] > 0.0) ||
                        ((phi[ii][ix][iy] == 0.0) && (phi[ii][ixp][iy] > 0.0) || (phi[ii][ixm][iy] > 0.0) || (phi[ii][ix][iyp] > 0.0) || (phi[ii][ix][iym] > 0.0)))
                    {
                        phinum++;
                        phiIdx[phinum][ix][iy] = ii;
                    }
                }
                phiNum[ix][iy] = phinum;
            }
        }
#pragma omp barrier
        // --------------------- Calculate concentration  --------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                if (ix == ndmx)
                {
                    ixp = 0;
                }
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (iy == ndmy)
                {
                    iyp = 0;
                }
                if (iy == 0)
                {
                    iym = ndmy;
                }
                if (phi[0][ix][iy] == 0.0)
                {
                    conp[1][ix][iy] = cont[ix][iy] * phi[1][ix][iy];
                    conp[2][ix][iy] = cont[ix][iy] * phi[2][ix][iy];
                    // Correct abnormal calculation at solid edge
                    if ((phi[0][ixp][iy] > 0.0) || (phi[0][ixm][iy] > 0.0) || (phi[0][ix][iyp] > 0.0) || (phi[0][ix][iym] > 0.0))
                    {
                        conp[1][ix][iy] = calC1e(temp[ix][iy]);
                        conp[2][ix][iy] = calC2e(temp[ix][iy]);
                    }
                    conp[0][ix][iy] = calC01e(temp[ix][iy]) * phi[1][ix][iy] + calC02e(temp[ix][iy]) * phi[2][ix][iy];
                }
                else if (phi[0][ix][iy] > 0.0 && phi[0][ix][iy] < 1.0)
                {
                    conp[1][ix][iy] = calC1e(temp[ix][iy]);
                    conp[2][ix][iy] = calC2e(temp[ix][iy]);
                    conp[0][ix][iy] = (cont[ix][iy] - conp[1][ix][iy] * phi[1][ix][iy] - conp[2][ix][iy] * phi[2][ix][iy]) / phi[0][ix][iy];
                    // Correct abnormal calculation at liquid edge
                    if (phi[0][ix][iy] < 0.05)
                    {
                        conp[0][ix][iy] = (calC01e(temp[ix][iy]) * phi[1][ix][iy] + calC02e(temp[ix][iy]) * phi[2][ix][iy]) / (phi[1][ix][iy] + phi[2][ix][iy]);
                    }
                    if (conp[0][ix][iy] > 1.0)
                    {
                        conp[0][ix][iy] = 1.0;
                    }
                    if (conp[0][ix][iy] < 0.0)
                    {
                        conp[0][ix][iy] = 0.0;
                    }
                }
                else if (phi[0][ix][iy] == 1.0)
                {
                    conp[1][ix][iy] = calC1e(temp[ix][iy]);
                    conp[2][ix][iy] = calC2e(temp[ix][iy]);
                    conp[0][ix][iy] = cont[ix][iy];
                }
                cont[ix][iy] = conp[1][ix][iy] * phi[1][ix][iy] + conp[2][ix][iy] * phi[2][ix][iy] + conp[0][ix][iy] * phi[0][ix][iy];
                if (cont[ix][iy] > 1.0)
                {
                    cont[ix][iy] = 1.0;
                }
                if (cont[ix][iy] < 0.0)
                {
                    cont[ix][iy] = 0.0;
                }
            }
        }
        // --------------------- Correct concentration in liquid phase for mass conservation --------------------------
#pragma omp barrier
        // collect mass
        if (th_id == 0)
        {
            sum0 = 0.0;
            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    sum0 += cont[ix][iy];
                }
            }
            c00 = sum0 / NDX / NDY;
        }
#pragma omp barrier
        // ---------------------------------  Evolution Equation of Concentration field ------------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                if (ix == ndmx)
                {
                    ixp = 0;
                }
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (iy == ndmy)
                {
                    iyp = 0;
                }
                if (iy == 0)
                {
                    iym = ndmy;
                }
                //拡散方程式内における微分計算
                for (ii = 1; ii < nm; ii++)
                {
                    sumcs = 0.25 * ((phi[ii][ixp][iy] - phi[ii][ixm][iy]) * (conp[ii][ixp][iy] - conp[ii][ixm][iy]) + (phi[ii][ix][iyp] - phi[ii][ix][iym]) * (conp[ii][ix][iyp] - conp[ii][ix][iym])) / dx / dx + phi[ii][ix][iy] * (conp[ii][ixp][iy] + conp[ii][ixm][iy] + conp[ii][ix][iyp] + conp[ii][ix][iym] - 4.0 * conp[ii][ix][iy]) / dx / dx;
                }
                sumcl = 0.25 * ((phi[0][ixp][iy] - phi[0][ixm][iy]) * (conp[0][ixp][iy] - conp[0][ixm][iy]) + (phi[0][ix][iyp] - phi[0][ix][iym]) * (conp[0][ix][iyp] - conp[0][ix][iym])) / dx / dx + phi[0][ix][iy] * (conp[0][ixp][iy] + conp[0][ixm][iy] + conp[0][ix][iyp] + conp[0][ix][iym] - 4.0 * conp[0][ix][iy]) / dx / dx;
                cddtt = Ds * sumcs + Dl * sumcl;
                cont2[ix][iy] = cont[ix][iy] + cddtt * dtime;
                // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
            }
        }

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                cont[ix][iy] = cont2[ix][iy];
            }
        }

        // cooling down
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                temp[ix][iy] -= dts;
            }
        }
        istep = istep + 1;
#pragma omp barrier
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
    int i, j;

    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/phi/2d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            fprintf(stream, "%e   ", phi[1][i][j]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream); //ファイルをクローズ

    FILE *streamc; //ストリームのポインタ設定
    char bufferc[30];
    sprintf(bufferc, "data/con/2d%d.csv", step);
    streamc = fopen(bufferc, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            fprintf(streamc, "%e   ", cont[i][j]);
            fprintf(streamc, "\n");
        }
    }
    fclose(streamc); //ファイルをクローズ

    FILE *streamt; //ストリームのポインタ設定
    char buffert[30];
    sprintf(buffert, "data/temp/2d%d.csv", step);
    streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            fprintf(streamt, "%e   ", temp[i][j] / 10.0);
            fprintf(streamt, "\n");
        }
        fclose(streamt); //ファイルをクローズ
    }
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
