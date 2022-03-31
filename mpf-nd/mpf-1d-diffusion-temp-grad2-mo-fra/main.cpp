
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <omp.h>

using namespace std;

#define N 2
#define NTH 4
#define ND 128
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;
int rows = ND / NTH;

int nstep = 200001;
int pstep = 2000;

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

double temp0 = 2.00;
double gradT = 0.00;
double dts = 0.00 / nstep;
double cl = 0.2;

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];
double phi[N][ND], phi2[N][ND];
double cont[ND], cont2[ND], conp[N][ND];
double temp[ND];
int phiNum[ND];
int phiIdx[N + 1][ND];

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

    for (i = 0; i <= ndm; i++)
    {
        temp[i] = temp0 + gradT * i;
    }

    sum0 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        if (i <= ND / 8)
        {
            phi[1][i] = 1.0;
            conp[1][i] = calC1e(temp[i]);
            phi[0][i] = 0.0;
            conp[0][i] = calC01e(temp[i]);
        }
        else
        {
            phi[1][i] = 0.0;
            conp[1][i] = calC1e(temp[i]);
            phi[0][i] = 1.0;
            conp[0][i] = cl;
        }
        cont[i] = conp[1][i] * phi[1][i] + conp[0][i] * phi[0][i];
        sum0 += cont[i];
    }
    c0 = sum0 / ND;
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

        phinum0 = 0;
        for (ni = 0; ni <= nm; ni++)
        {
            if ((phi[ni][i] > 0.0) ||
                ((phi[ni][i] == 0.0) && (phi[ni][ip] > 0.0) || (phi[ni][im] > 0.0)))
            {
                phinum0++;
                phiIdx[phinum0][i] = ni;
            }
        }
        phiNum[i] = phinum0;
    }

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

        if ((istep % pstep == 0) && (th_id == 0))
        {
            datasave(istep);
            cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
            cout << "The nominal concnetration is " << c00 << endl;
            cout << "The th_id is " << th_id << endl;
        }

        // ---------------------------------  Evolution Equation of Phase fields ------------------------------------
        for (ix = start; ix <= end; ix++)
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
                        dF = calDF10(conp[0][ix], temp[ix], S0);
                    }
                    else if (ii == 0 && jj == 1)
                    {
                        dF = -calDF10(conp[0][ix], temp[ix], S0);
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

        for (ix = start; ix <= end; ix++)
        {
            for (kk = 0; kk <= nm; kk++)
            {
                phi[kk][ix] = phi2[kk][ix];
            }
        }

        //
        for (ix = start; ix <= end; ix++)
        {
            psum = 0.0;
            for (kk = 0; kk <= nm; kk++)
            {
                psum += phi[kk][ix];
            }
            for (kk = 0; kk <= nm; kk++)
            {
                phi[kk][ix] = phi[kk][ix] / psum;
                phi2[kk][ix] = phi[kk][ix];
            }
        }
#pragma omp barrier
        // ---------------------------  Collect information of phase fields ----------------------------------
        for (ix = start; ix <= end; ix++)
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
#pragma omp barrier
        // --------------------- Calculate concentration in the interface and liquid phase --------------------------
        for (ix = start; ix <= end; ix++)
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
            if (phi[0][ix] == 0.0)
            {
                conp[1][ix] = cont[ix];
                // Correct abnormal calculation at solid edge
                if ((phi[0][ixp] > 0.0) || (phi[0][ixm] > 0.0))
                {
                    conp[1][ix] = calC1e(temp[ix]);
                }
                conp[0][ix] = calC01e(temp[ix]);
            }
            else if (phi[0][ix] > 0.0 && phi[0][ix] < 1.0)
            {
                conp[1][ix] = calC1e(temp[ix]);
                conp[0][ix] = (cont[ix] - conp[1][ix] * phi[1][ix]) / phi[0][ix];
                // Correct abnormal calculation at liquid edge
                if (phi[0][ix] < 0.05)
                {
                    conp[0][ix] = calC01e(temp[ix]);
                }
                if (conp[0][ix] > 1.0)
                {
                    conp[0][ix] = 1.0;
                }
                if (conp[0][ix] < 0.0)
                {
                    conp[0][ix] = 0.0;
                }
            }
            else if (phi[0][ix] == 1.0)
            {
                conp[1][ix] = calC1e(temp[ix]);
                conp[0][ix] = cont[ix];
            }
            cont[ix] = conp[1][ix] * phi[1][ix] + conp[0][ix] * phi[0][ix];
            if (cont[ix] > 1.0)
            {
                cont[ix] = 1.0;
            }
            if (cont[ix] < 0.0)
            {
                cont[ix] = 0.0;
            }
        }
        // --------------------- Correct concentration in liquid phase for mass conservation --------------------------
#pragma omp barrier
        // collect mass
        if (th_id == 0)
        {
            sum0 = 0.0;
            for (ix = 0; ix <= ndm; ix++)
            {
                sum0 += cont[ix];
            }
            c00 = sum0 / ND;
            dc0 = sum0 / ND - c0;
            // correction for mass conservation in the interface region
            for (ix = 0; ix <= ndm; ix++)
            {
                cont[ix] = cont[ix] - dc0;
                if (cont[ix] > 1.0)
                {
                    cont[ix] = 1.0;
                }
                if (cont[ix] < 0.0)
                {
                    cont[ix] = 0.0;
                }
            }
        }
#pragma omp barrier
        // ---------------------------------  Evolution Equation of Concentration field ------------------------------------
        for (ix = start; ix <= end; ix++)
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
            dev1_s = 0.25 * ((phi[1][ixp] - phi[1][ixm]) * (conp[1][ixp] - conp[1][ixm])) / dx / dx;
            dev1_l = 0.25 * ((phi[0][ixp] - phi[0][ixm]) * (conp[0][ixp] - conp[0][ixm])) / dx / dx;
            dev2_s = phi[1][ix] * (conp[1][ixp] + conp[1][ixm] - 2.0 * conp[1][ix]) / dx / dx;
            dev2_l = phi[0][ix] * (conp[0][ixp] + conp[0][ixm] - 2.0 * conp[0][ix]) / dx / dx;

            cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l); //拡散方程式[式(4.42)]
            cont2[ix] = cont[ix] + cddtt * dtime;                    //濃度場の時間発展(陽解法)
                                                                     // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
        }

        for (ix = start; ix <= end; ix++)
        {
            cont[ix] = cont2[ix];
        }

        // cooling down
        for (ix = start; ix <= end; ix++)
        {
            temp[ix] -= dts;
        }
        // --------------------------------  Long-range solute diffusion calculation -----------------------------------

        // ----------------------------------------------  Moving frame  -----------------------------------------------
#pragma omp barrier
        if (th_id == 0)
        {
            if (istep % pstep == 0)
            {
                cout << "the thread which is moving frame is " << th_id << endl;
            }
            intpos = 0;
            // check if the bottom is solid
            if (phi[0][0] != 0.0)
            {
                hasS = 0;
            }
            else if (phi[0][0] == 0.0)
            {
                hasS = 1;
            }
            // search interface front
            if (hasS == 1)
            {
                allS = 1;
                for (ix = 0; ix <= ndm; ix++)
                {
                    if (phi[0][ix] > 0.0)
                    {
                        allS = 0;
                    }
                    if (allS == 0)
                    {
                        intpos = ix;
                        break;
                    }
                }

                allL = 0;
                for (ix = intpos; ix <= ndm; ix++)
                {
                    if (phi[0][ix] == 0.0)
                    {
                        allL = 1;
                    }
                    if (allL == 1)
                    {
                        intpos = ix;
                        break;
                    }
                }
            }
            // check the distance from the middle of the domain
            if (intpos > ND / 4)
            {
                dist = intpos - ND / 4;
                cout << "The dist is " << dist << endl;
                for (ix = 0; ix <= (ndm - dist); ix++)
                {
                    // temp
                    temp[ix] = temp[ix + dist];
                    // cont
                    cont[ix] = cont[ix + dist];
                    // phi
                    phi[0][ix] = phi[0][ix + dist];
                    phi[1][ix] = phi[1][ix + dist];
                    // conp
                    conp[0][ix] = conp[0][ix + dist];
                    conp[1][ix] = conp[1][ix + dist];
                }
                for (ix = (ndm - dist + 1); ix <= ndm; ix++)
                {
                    // temp
                    temp[ix] = temp[ndm - dist] + gradT * (ix - ndm + dist);
                    // cont
                    cont[ix] = cont[ndm - dist];
                    // phi
                    phi[0][ix] = 1.0;
                    phi[1][ix] = 0.0;
                    // conp
                    conp[0][ix] = cont[ix];
                    conp[1][ix] = calC1e(temp[ix]);
                }
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
    int i;

    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/phi/1d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndm; i++)
    {
        fprintf(stream, "%e   ", phi[1][i]);
        fprintf(stream, "\n");
    }
    fclose(stream); //ファイルをクローズ

    FILE *streamc; //ストリームのポインタ設定
    char bufferc[30];
    sprintf(bufferc, "data/con/1d%d.csv", step);
    streamc = fopen(bufferc, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndm; i++)
    {
        fprintf(streamc, "%e   ", cont[i]);
        fprintf(streamc, "\n");
    }
    fclose(streamc); //ファイルをクローズ

    FILE *streamcl; //ストリームのポインタ設定
    char buffercl[30];
    sprintf(buffercl, "data/conl/1d%d.csv", step);
    streamcl = fopen(buffercl, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndm; i++)
    {
        fprintf(streamcl, "%e   ", conp[0][i]);
        fprintf(streamcl, "\n");
    }
    fclose(streamcl); //ファイルをクローズ

    FILE *streamcs; //ストリームのポインタ設定
    char buffercs[30];
    sprintf(buffercs, "data/cons/1d%d.csv", step);
    streamcs = fopen(buffercs, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndm; i++)
    {
        fprintf(streamcs, "%e   ", conp[1][i]);
        fprintf(streamcs, "\n");
    }
    fclose(streamcs); //ファイルをクローズ

    FILE *streamt; //ストリームのポインタ設定
    char buffert[30];
    sprintf(buffert, "data/temp/1d%d.csv", step);
    streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndm; i++)
    {
        fprintf(streamt, "%e   ", temp[i] / 10.0);
        fprintf(streamt, "\n");
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
