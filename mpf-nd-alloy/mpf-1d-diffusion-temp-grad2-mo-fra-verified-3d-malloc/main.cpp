#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <omp.h>
#include "CImg.h" //CImg ライブラリ（描画用）使用のためのヘッダ

using namespace std;
using namespace cimg_library;

#define N 3
#define NTH 8
#define NDX 64
#define NDY 64
#define NDZ 64
#define NDL 2560
#define PI 3.14159

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int ndml = NDL - 1;
int mid = NDX / 2;
int rows = NDX / NTH;
int rowsl = NDL / NTH;

int nstep = 2001;
int pstep = 2000;

double dx = 1.0;
double dtime = 1.0;

double gamma0 = 0.1;
double mobi = 0.25;
double delta = 5.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double S0 = 0.03;

double Dl = 0.1;
double Ds = 2.0e-4;

double gradT = 0.002;
// double rateT = 0.000002;
// double temp0 = -0.60 - NDZ / 4 * gradT;
double rateT = 0.000006;
double temp0 = -1.30 - NDZ / 4 * gradT;
double cl = 0.5;

double alpha_d = dtime * Dl / dx / dx;
double alpha_m = dtime / dx / dx * mobi * A0;

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];
// double phi[N][NDX][NDY][NDZ], phi2[N][NDX][NDY][NDZ];
// double cont[NDX][NDY][NDZ], cont2[NDX][NDY][NDZ], conp[N][NDX][NDY][NDZ];
double conlr[NDL], conlr2[NDL];
// double (*temp)[NDX][NDY][NDZ];
// int phiNum[NDX][NDY][NDZ];
// int phiIdx[N + 1][NDX][NDY][NDZ];

CImg<unsigned char> ch_fld(NDX, NDZ, 1, 3);
char outFileCh_xz[64];
char outFiles_xz[64];

void datasave(int step);
double calC01e(double temp0), calC1e(double temp0), calC02e(double temp0), calC2e(double temp0);
double calDF10(double con0, double temp0, double dS), calDF20(double con0, double temp0, double dS);

int main(void)
{

    int i, j, k, im, ip, jm, jp, km, kp;
    int ni, phinum0;
    int intpos, dist, hasS, allS, allL;
    double c0, c00, dc0, sum0, sumplane;

    double(*phi)[N][NDX][NDY][NDZ] = malloc(sizeof(*phi));
    double(*phi2)[N][NDX][NDY][NDZ] = malloc(sizeof(*phi2));
    int(*phiIdx)[N + 1][NDX][NDY][NDZ] = malloc(sizeof(*phiIdx));
    int(*phiNum)[NDX][NDY][NDZ] = malloc(sizeof(*phiNum));
    double(*cont)[NDX][NDY][NDZ] = malloc(sizeof(*cont));
    double(*cont2)[NDX][NDY][NDZ] = malloc(sizeof(*cont2));
    double(*conp)[N][NDX][NDY][NDZ] = malloc(sizeof(*conp));
    double(*temp)[NDX][NDY][NDZ] = malloc(sizeof(*temp));

    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "concenration field stablity number is: " << alpha_d << endl;
    cout << "phase field stability number is: " << alpha_m << endl;

    if ((alpha_d > 0.15) || (alpha_m > 0.15))
    {
        cout << "The computation is unstable, please change input parameters!" << endl;
        goto terminal;
    }

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

    mij[1][2] = M0 * 0.1;
    mij[2][1] = M0 * 0.1;

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                (*temp)[i][j][k] = temp0 + gradT * i * dx;
            }
        }
    }

    sum0 = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {

                // if ((i < NDX / 2) && (k < NDZ / 4))
                if (((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < (NDX * NDX / 2.0 / PI)) && (k < NDZ / 4))
                {
                    (*phi)[1][i][j][k] = 1.0;
                    (*conp)[1][i][j][k] = calC1e((*temp)[i][j][k]);
                    (*phi)[2][i][j][k] = 0.0;
                    (*conp)[2][i][j][k] = calC2e((*temp)[i][j][k]);
                    (*phi)[0][i][j][k] = 0.0;
                    (*conp)[0][i][j][k] = calC01e((*temp)[i][j][k]);
                }
                else if (((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) >= (NDX * NDX / 2.0 / PI)) && (k < NDZ / 4))
                {
                    (*phi)[1][i][j][k] = 0.0;
                    (*conp)[1][i][j][k] = calC1e((*temp)[i][j][k]);
                    (*phi)[2][i][j][k] = 1.0;
                    (*conp)[2][i][j][k] = calC2e((*temp)[i][j][k]);
                    (*phi)[0][i][j][k] = 0.0;
                    (*conp)[0][i][j][k] = calC02e((*temp)[i][j][k]);
                }
                else
                {
                    (*phi)[1][i][j][k] = 0.0;
                    (*conp)[1][i][j][k] = calC1e((*temp)[i][j][k]);
                    (*phi)[2][i][j][k] = 0.0;
                    (*conp)[2][i][j][k] = calC2e((*temp)[i][j][k]);
                    (*phi)[0][i][j][k] = 1.0;
                    (*conp)[0][i][j][k] = cl;
                }
                (*cont)[i][j][k] = (*conp)[1][i][j][k] * (*phi)[1][i][j][k] + (*conp)[2][i][j][k] * (*phi)[2][i][j][k] + (*conp)[0][i][j][k] * (*phi)[0][i][j][k];
                sum0 += (*cont)[i][j][k];
            }
        }
    }
    c0 = sum0 / NDX / NDY / NDZ;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                if (i == ndmx)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndmx;
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                if (k == ndmz)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndmz;
                }

                phinum0 = 0;
                for (ni = 0; ni <= nm; ni++)
                {
                    if (((*phi)[ni][i][j][k] > 0.0) ||
                        (((*phi)[ni][i][j][k] == 0.0) && ((*phi)[ni][ip][j][k] > 0.0) ||
                         ((*phi)[ni][im][j][k] > 0.0) ||
                         ((*phi)[ni][i][jp][k] > 0.0) ||
                         ((*phi)[ni][i][jm][k] > 0.0) ||
                         ((*phi)[ni][i][j][kp] > 0.0) ||
                         ((*phi)[ni][i][j][km] > 0.0)))
                    {
                        phinum0++;
                        (*phiIdx)[phinum0][i][j][k] = ni;
                    }
                }
                (*phiNum)[i][j][k] = phinum0;
            }
        }
    }
#pragma omp parallel num_threads(NTH)
    {
        int istep, th_id;
        int start, end, offset;
        int startl, endL, offsetl;
        int ix, ixm, ixp, iy, iym, iyp, iz, izm, izp;
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
            cout << "--" << endl;
            cout << "   The nominal concnetration is " << c00 << endl;
            // ****** YZ *******
            cimg_forXY(ch_fld, x, z)
            {
                ch_fld(x, z, 0) = 255. * ((*cont)[x][NDY / 2][z]); // red
                ch_fld(x, z, 1) = 255. * ((*cont)[x][NDY / 2][z]); // green
                ch_fld(x, z, 2) = 255. * ((*cont)[x][NDY / 2][z]); // blue
            }
            sprintf(outFileCh_xz, "figures/con/2d%d.png", istep); // generate imagefile
            ch_fld.save_jpeg(outFileCh_xz);                       // save imagegilee
        }

        // ---------------------------------  Evolution Equation of Phase fields ------------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
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
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }

                    for (n1 = 1; n1 <= (*phiNum)[ix][iy][iz]; n1++)
                    {
                        ii = (*phiIdx)[n1][ix][iy][iz];
                        pddtt = 0.0;
                        for (n2 = 1; n2 <= (*phiNum)[ix][iy][iz]; n2++)
                        {
                            jj = (*phiIdx)[n2][ix][iy][iz];
                            dsum = 0.0;
                            for (n3 = 1; n3 <= (*phiNum)[ix][iy][iz]; n3++)
                            {
                                kk = (*phiIdx)[n3][ix][iy][iz];

                                termiikk = aij[ii][kk] * ((*phi)[kk][ixp][iy][iz] + (*phi)[kk][ixm][iy][iz] + (*phi)[kk][ix][iyp][iz] + (*phi)[kk][ix][iym][iz] + (*phi)[kk][ix][iy][izp] + (*phi)[kk][ix][iy][izm] - 6.0 * (*phi)[kk][ix][iy][iz]) / (dx * dx);

                                termjjkk = aij[jj][kk] * ((*phi)[kk][ixp][iy][iz] + (*phi)[kk][ixm][iy][iz] + (*phi)[kk][ix][iyp][iz] + (*phi)[kk][ix][iym][iz] + (*phi)[kk][ix][iy][izp] + (*phi)[kk][ix][iy][izm] - 6.0 * (*phi)[kk][ix][iy][iz]) / (dx * dx);

                                dsum += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * (*phi)[kk][ix][iy][iz];
                            }
                            if (ii == 1 && jj == 0)
                            {
                                dF = calDF10((*conp)[0][ix][iy][iz], (*temp)[ix][iy][iz], S0);
                            }
                            else if (ii == 0 && jj == 1)
                            {
                                dF = -calDF10((*conp)[0][ix][iy][iz], (*temp)[ix][iy][iz], S0);
                            }
                            else if (ii == 2 && jj == 0)
                            {
                                dF = calDF20((*conp)[0][ix][iy][iz], (*temp)[ix][iy][iz], S0);
                            }
                            else if (ii == 0 && jj == 2)
                            {
                                dF = -calDF20((*conp)[0][ix][iy][iz], (*temp)[ix][iy][iz], S0);
                            }
                            else
                            {
                                dF = 0.0;
                            }
                            pddtt += -2.0 * mij[ii][jj] / double((*phiNum)[ix][iy][iz]) * (dsum - 8.0 / PI * dF * sqrt((*phi)[ii][ix][iy][iz] * (*phi)[jj][ix][iy][iz]));
                        }
                        (*phi2)[ii][ix][iy][iz] = (*phi)[ii][ix][iy][iz] + pddtt * dtime;
                        if ((*phi2)[ii][ix][iy][iz] >= 1.0)
                        {
                            (*phi2)[ii][ix][iy][iz] = 1.0;
                        }
                        if ((*phi2)[ii][ix][iy][iz] <= 0.0)
                        {
                            (*phi2)[ii][ix][iy][iz] = 0.0;
                        }
                    }
                } // ix
            }     // iy
        }         // iz

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    for (kk = 0; kk <= nm; kk++)
                    {
                        (*phi)[kk][ix][iy][iz] = (*phi2)[kk][ix][iy][iz];
                    }
                }
            }
        }

        //
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    psum = 0.0;
                    for (kk = 0; kk <= nm; kk++)
                    {
                        psum += (*phi)[kk][ix][iy][iz];
                    }
                    for (kk = 0; kk <= nm; kk++)
                    {
                        (*phi)[kk][ix][iy][iz] = (*phi)[kk][ix][iy][iz] / psum;
                        (*phi2)[kk][ix][iy][iz] = (*phi)[kk][ix][iy][iz];
                    }
                }
            }
        }
#pragma omp barrier
        // ---------------------------  Collect information of phase fields ----------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
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
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }

                    phinum = 0;
                    for (ii = 0; ii <= nm; ii++)
                    {
                        if (((*phi)[ii][ix][iy][iz] > 0.0) ||
                            (((*phi)[ii][ix][iy][iz] == 0.0) &&
                                 ((*phi)[ii][ixp][iy][iz] > 0.0) ||
                             ((*phi)[ii][ixm][iy][iz] > 0.0) ||
                             ((*phi)[ii][ix][iyp][iz] > 0.0) ||
                             ((*phi)[ii][ix][iym][iz] > 0.0) ||
                             ((*phi)[ii][ix][iy][izp] > 0.0) ||
                             ((*phi)[ii][ix][iy][izm] > 0.0)))
                        {
                            phinum++;
                            (*phiIdx)[phinum][ix][iy][iz] = ii;
                        }
                    }
                    (*phiNum)[ix][iy][iz] = phinum;
                }
            }
        }
#pragma omp barrier
        // --------------------- Calculate concentration  --------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
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
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }
                    if ((*phi)[0][ix][iy][iz] == 0.0)
                    {
                        if ((*phi)[1][ix][iy][iz] == 1.0)
                        {
                            (*conp)[1][ix][iy][iz] = (*cont)[ix][iy][iz];
                        }
                        else
                        {
                            (*conp)[1][ix][iy][iz] = calC1e((*temp)[ix][iy][iz]);
                        }
                        if ((*phi)[2][ix][iy][iz] == 1.0)
                        {
                            (*conp)[2][ix][iy][iz] = (*cont)[ix][iy][iz];
                        }
                        else
                        {
                            (*conp)[2][ix][iy][iz] = calC2e((*temp)[ix][iy][iz]);
                        }
                        // Correct abnormal calculation at solid edge
                        if (((*phi)[0][ixp][iy][iz] > 0.0) || ((*phi)[0][ixm][iy][iz] > 0.0) || ((*phi)[0][ix][iyp][iz] > 0.0) || ((*phi)[0][ix][iym][iz] > 0.0) || ((*phi)[0][ix][iy][izp] > 0.0) || ((*phi)[0][ix][iy][izm] > 0.0))
                        {
                            (*conp)[1][ix][iy][iz] = calC1e((*temp)[ix][iy][iz]);
                            (*conp)[2][ix][iy][iz] = calC2e((*temp)[ix][iy][iz]);
                        }
                        (*conp)[0][ix][iy][iz] = calC01e((*temp)[ix][iy][iz]) * (*phi)[1][ix][iy][iz] + calC02e((*temp)[ix][iy][iz]) * (*phi)[2][ix][iy][iz];
                    }
                    else if ((*phi)[0][ix][iy][iz] > 0.0 && (*phi)[0][ix][iy][iz] < 1.0)
                    {
                        (*conp)[1][ix][iy][iz] = calC1e((*temp)[ix][iy][iz]);
                        (*conp)[2][ix][iy][iz] = calC2e((*temp)[ix][iy][iz]);
                        (*conp)[0][ix][iy][iz] = ((*cont)[ix][iy][iz] - (*conp)[1][ix][iy][iz] * (*phi)[1][ix][iy][iz] - (*conp)[2][ix][iy][iz] * (*phi)[2][ix][iy][iz]) / (*phi)[0][ix][iy][iz];
                        // Correct abnormal calculation at liquid edge
                        if ((*phi)[0][ix][iy][iz] < 0.05)
                        {
                            (*conp)[0][ix][iy][iz] = (calC01e((*temp)[ix][iy][iz]) * (*phi)[1][ix][iy][iz] + calC02e((*temp)[ix][iy][iz]) * (*phi)[2][ix][iy][iz]) / ((*phi)[1][ix][iy][iz] + (*phi)[2][ix][iy][iz]);
                        }
                        if ((*conp)[0][ix][iy][iz] > 1.0)
                        {
                            (*conp)[0][ix][iy][iz] = 1.0;
                        }
                        if ((*conp)[0][ix][iy][iz] < 0.0)
                        {
                            (*conp)[0][ix][iy][iz] = 0.0;
                        }
                    }
                    else if ((*phi)[0][ix][iy][iz] == 1.0)
                    {
                        (*conp)[1][ix][iy][iz] = calC1e((*temp)[ix][iy][iz]);
                        (*conp)[2][ix][iy][iz] = calC2e((*temp)[ix][iy][iz]);
                        (*conp)[0][ix][iy][iz] = (*cont)[ix][iy][iz];
                    }
                    (*cont)[ix][iy][iz] = (*conp)[1][ix][iy][iz] * (*phi)[1][ix][iy][iz] + (*conp)[2][ix][iy][iz] * (*phi)[2][ix][iy][iz] + (*conp)[0][ix][iy][iz] * (*phi)[0][ix][iy][iz];
                    if ((*cont)[ix][iy][iz] > 1.0)
                    {
                        (*cont)[ix][iy][iz] = 1.0;
                    }
                    if ((*cont)[ix][iy][iz] < 0.0)
                    {
                        (*cont)[ix][iy][iz] = 0.0;
                    }
                }
            }
        }
        // --------------------- Correct concentration in liquid phase for mass conservation --------------------------
#pragma omp barrier
        // collect mass
        if (th_id == 0 && (istep > nstep / 5))
        {
            sum0 = 0.0;
            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    for (iz = 0; iz <= ndmz; iz++)
                    {
                        sum0 += (*cont)[ix][iy][iz];
                    }
                }
            }
            c00 = sum0 / NDX / NDY / NDZ;
            dc0 = c00 - cl;
            if (dc0 >= 0.003)
            {
                for (ix = 0; ix <= ndmx; ix++)
                {
                    for (iy = 0; iy <= ndmy; iy++)
                    {
                        for (iz = 0; iz <= ndmz; iz++)
                        {
                            (*cont)[ix][iy][iz] = (*cont)[ix][iy][iz] - dc0;
                            if ((*cont)[ix][iy][iz] > 1.0)
                            {
                                (*cont)[ix][iy][iz] = 1.0;
                            }
                            if ((*cont)[ix][iy][iz] < 0.0)
                            {
                                (*cont)[ix][iy][iz] = 0.0;
                            }
                        }
                    }
                }
            }
        }
#pragma omp barrier
        // ---------------------------------  Evolution Equation of Concentration field ------------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
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
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }
                    //拡散方程式内における微分計算
                    for (ii = 1; ii < nm; ii++)
                    {
                        sumcs = 0.25 * (((*phi)[ii][ixp][iy][iz] - (*phi)[ii][ixm][iy][iz]) * ((*conp)[ii][ixp][iy][iz] - (*conp)[ii][ixm][iy][iz]) + ((*phi)[ii][ix][iyp][iz] - (*phi)[ii][ix][iym][iz]) * ((*conp)[ii][ix][iyp][iz] - (*conp)[ii][ix][iym][iz]) + ((*phi)[ii][ix][iy][izp] - (*phi)[ii][ix][iy][izm]) * ((*conp)[ii][ix][iy][izp] - (*conp)[ii][ix][iy][izm])) / dx / dx +
                                (*phi)[ii][ix][iy][iz] * ((*conp)[ii][ixp][iy][iz] + (*conp)[ii][ixm][iy][iz] + (*conp)[ii][ix][iyp][iz] + (*conp)[ii][ix][iym][iz] + (*conp)[ii][ix][iy][izp] + (*conp)[ii][ix][iy][izm] - 6.0 * (*conp)[ii][ix][iy][iz]) / dx / dx;
                    }
                    sumcl = 0.25 * (((*phi)[0][ixp][iy][iz] - (*phi)[0][ixm][iy][iz]) * ((*conp)[0][ixp][iy][iz] - (*conp)[0][ixm][iy][iz]) + ((*phi)[0][ix][iyp][iz] - (*phi)[0][ix][iym][iz]) * ((*conp)[0][ix][iyp][iz] - (*conp)[0][ix][iym][iz]) + ((*phi)[0][ix][iy][izp] - (*phi)[0][ix][iy][izm]) * ((*conp)[0][ix][iy][izp] - (*conp)[0][ix][iy][izm])) / dx / dx +
                            (*phi)[0][ix][iy][iz] * ((*conp)[0][ixp][iy][iz] + (*conp)[0][ixm][iy][iz] + (*conp)[0][ix][iyp][iz] + (*conp)[0][ix][iym][iz] + (*conp)[0][ix][iy][izp] + (*conp)[0][ix][iy][izm] - 6.0 * (*conp)[0][ix][iy][iz]) / dx / dx;
                    cddtt = Ds * sumcs + Dl * sumcl;
                    (*cont2)[ix][iy][iz] = (*cont)[ix][iy][iz] + cddtt * dtime;
                    // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
                }
            }
        }

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    (*cont)[ix][iy][iz] = (*cont2)[ix][iy][iz];
                }
            }
        }

        // cooling down
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    (*temp)[ix][iy][iz] -= rateT * dtime;
                }
            }
        }

        //----------------------------------------------  Moving frame  -----------------------------------------------
#pragma omp barrier
        if (th_id == 0)
        {
            // check if the bottom is solid
            sumplane = 0.0;
            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    if ((*phi)[0][ix][iy][0] != 0.0)
                    {
                        hasS = 0;
                    }
                    sumplane += (*phi)[0][ix][iy][0];
                    if ((sumplane == 0.0) && (iy == ndmy) && (ix == ndmx))
                    {
                        hasS = 1;
                    }
                }
            }
            // search interface front
            intpos = 0;
            if (hasS == 1)
            {
                allS = 1;
                for (iz = 0; iz <= ndmz; iz++)
                {
                    if (allS == 0)
                    {
                        intpos = iz - 1;
                        break;
                    }
                    for (ix = 0; ix <= ndmx; ix++)
                    {
                        for (iy = 0; iy <= ndmy; iy++)
                        {
                            if ((*phi)[0][ix][iy][iz] > 0.0)
                            {
                                allS = 0;
                                break;
                            }
                        }
                        if (allS == 0)
                        {
                            break;
                        }
                    }
                }

                allL = 0;
                for (iz = intpos; iz <= ndmx; iz++)
                {
                    sumplane = 0.0;
                    for (ix = 0; ix <= ndmx; ix++)
                    {
                        for (iy = 0; iy <= ndmy; iy++)
                        {
                            sumplane += (*phi)[0][ix][iy][iz];
                        }
                    }
                    if (sumplane == double(NDX * NDY))
                    {
                        allL = 1;
                    }
                    if (allL == 1)
                    {
                        intpos = iz;
                        break;
                    }
                }
            }
            // check the distance from the middle of the domain
            if (intpos > mid)
            {
                dist = intpos - mid;
                cout << "--" << endl;
                cout << "    the distance away from middle is " << dist << endl;
                cout << "--" << endl;
                cout << "    the interface temperature is " << (*temp)[NDX / 2][NDY / 2][mid] << endl;
                for (iz = 0; iz <= (ndmz - dist); iz++)
                {
                    for (ix = 0; ix <= ndmx; ix++)
                    {
                        for (iy = 0; iy <= ndmy; iy++)
                        {
                            // temp
                            (*temp)[ix][iy][iz] = (*temp)[ix][iy][iz + dist];
                            // cont
                            (*cont)[ix][iy][iz] = (*cont)[ix][iy][iz + dist];
                            // phi
                            (*phi)[0][ix][iy][iz] = (*phi)[0][ix][iy][iz + dist];
                            (*phi)[1][ix][iy][iz] = (*phi)[1][ix][iy][iz + dist];
                            (*phi)[2][ix][iy][iz] = (*phi)[2][ix][iy][iz + dist];
                            // conp
                            (*conp)[0][ix][iy][iz] = (*conp)[0][ix][iy][iz + dist];
                            (*conp)[1][ix][iy][iz] = (*conp)[1][ix][iy][iz + dist];
                            (*conp)[2][ix][iy][iz] = (*conp)[2][ix][iy][iz + dist];
                        }
                    }
                }
                for (iz = (ndmz - dist + 1); iz <= ndmz; iz++)
                {
                    for (ix = 0; ix <= ndmx; ix++)
                    {
                        for (iy = 0; iy <= ndmy; iy++)
                        {
                            // temp
                            (*temp)[ix][iy][iz] = (*temp)[ix][iy][ndmz - dist] + gradT * (iz - ndmz + dist) * dx;
                            // cont
                            (*cont)[ix][iy][iz] = cl;
                            // phi
                            (*phi)[0][ix][iy][iz] = 1.0;
                            (*phi)[1][ix][iy][iz] = 0.0;
                            (*phi)[2][ix][iy][iz] = 0.0;
                            // conp
                            (*conp)[0][ix][iy][iz] = cl;
                            (*conp)[1][ix][iy][iz] = calC1e((*temp)[ix][iy][iz]);
                            (*conp)[2][ix][iy][iz] = calC2e((*temp)[ix][iy][iz]);
                        }
                    }
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
terminal:;
    return 0;
}

void datasave(int step)
{
    int i, j, k;

    FILE *streamc;
    char bufferc[30];
    sprintf(bufferc, "data/con/3d%d.vtk", step);
    streamc = fopen(bufferc, "a");

    fprintf(streamc, "# vtk DataFile Version 1.0\n");
    fprintf(streamc, "phi_%d.vtk\n", step);
    fprintf(streamc, "ASCII\n");
    fprintf(streamc, "DATASET STRUCTURED_POINTS\n");
    fprintf(streamc, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
    fprintf(streamc, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(streamc, "ASPECT_RATIO 1.0 1.0 1.0\n");
    fprintf(streamc, "\n");
    fprintf(streamc, "POINT_DATA %d\n", NDX * NDY * NDZ);
    fprintf(streamc, "SCALARS scalars float\n");
    fprintf(streamc, "LOOKUP_TABLE default\n");

    for (k = 0; k <= ndmz; k++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (i = 0; i <= ndmx; i++)
            {
                fprintf(streamc, "%e\n", (*cont)[i][j][k]);
            }
        }
    }
    fclose(streamc);
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
