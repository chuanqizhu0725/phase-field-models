#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "CImg.h" //CImg ライブラリ（描画用）使用のためのヘッダ

using namespace std;
using namespace cimg_library;

#define PI 3.141592
#define RR 8.3145

void PhaseProperties(double delta, double gamma0, double mobi,
                     double A0, double W0, double M0,
                     double **aij, double **wij, double **mij,
                     double **anij, double **thij, double **vpij, double **etaij,
                     int ni, int nj, int nm)
{
    A0 = 8.0 * delta * gamma0 / PI / PI;
    W0 = 4.0 * gamma0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);

    for (ni = 0; ni <= nm; ni++)
    {
        for (nj = 0; nj <= nm; nj++)
        {
            wij[ni][nj] = W0;
            aij[ni][nj] = A0;
            mij[ni][nj] = M0;
            anij[ni][nj] = 0.0;
            thij[ni][nj] = 0.0;
            vpij[ni][nj] = 0.0;
            etaij[ni][nj] = 0.0;
            if (ni == nj)
            {
                wij[ni][nj] = 0.0;
                aij[ni][nj] = 0.0;
                mij[ni][nj] = 0.0;
            }
        }
    }
}

void TemperatureGradient(double ***temp, double temp0, double T_left, double Tg,
                         int ndmx, int ndmy, int ndmz, double dx,
                         int i, int j, int k)
{
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                temp[i][j][k] = T_left + i * dx * Tg;
            }
        }
    }
}

void CenterSeed(double ****phi,
                int NDX, int NDY, int NDZ, int ndmx, int ndmy, int ndmz,
                int i, int j, int k, int nm,
                double r0, double r, int xx0, int yy0, int zz0)
{
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                // r = sqrt((i - xx0) * (i - xx0) + (j - yy0) * (j - yy0) + (k - zz0) * (k - zz0));
                if (i <= NDX / 4)
                {
                    phi[nm][i][j][k] = 1.0;
                    phi[0][i][j][k] = 0.0;
                }
                else
                {
                    phi[0][i][j][k] = 1.0;
                    phi[nm][i][j][k] = 0.0;
                }
            }
        }
    }
}

void CollectPhaseFields(double ****phi, int ***phiNum, int ****phiIdx,
                        int start, int end, int ndmx, int ndmy, int ndmz,
                        int ix, int iy, int iz, int ixp, int ixm, int iyp, int iym, int izp, int izm,
                        int phinum, int nm)
{
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
                if (ix == 0)
                {
                    ixm = ndmx;
                }
                if (ix == ndmx)
                {
                    ixp = 0;
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
                    izp = 0;
                }
                if (iz == 0)
                {
                    izm = ndmz;
                }
                phinum = 0;
                for (int ii = 0; ii <= nm; ii++)
                {
                    if ((phi[ii][ix][iy][iz] > 0.0) ||
                        ((phi[ii][ix][iy][iz] == 0.0) && (phi[ii][ixp][iy][iz] > 0.0) ||
                         (phi[ii][ixm][iy][iz] > 0.0) ||
                         (phi[ii][ix][iyp][iz] > 0.0) ||
                         (phi[ii][ix][iym][iz] > 0.0) ||
                         (phi[ii][ix][iy][izp] > 0.0) ||
                         (phi[ii][ix][iy][izm] > 0.0)))
                    {
                        phinum++;
                        phiIdx[phinum][ix][iy][iz] = ii;
                    }
                }
                phiNum[ix][iy][iz] = phinum;
            }
        }
    }
}

void ComputePhaseFields(double ****phi, double ****phi2, int ***phiNum, int ****phiIdx,
                        double **aij, double **wij, double **mij, double **sij,
                        double **anij, double **thij, double **vpij, double **etaij, double astre, double astrem,
                        int start, int end, int ndmx, int ndmy, int ndmz, double dtime, double dx,
                        int ix, int iy, int iz, int ixp, int ixm, int iyp, int iym, int izp, int izm,
                        int ii, int jj, int kk, int n1, int n2, int n3, int nm,
                        double pddtt, double intsum, double psum,
                        double ****conp, double ***temp, double dH, double Tm, double sph_s)
{

    // *************  For anisotropy calculation **************
    double phidx, phidy, phidz;
    double phidxii, phidyii, phidzii;
    double phidxx, phidyy, phidzz;
    double phidxy, phidxz, phidyz;
    double phiabs, phiabsii;

    double th, vp, eta;
    double thii, vpii, etaii;
    double epsilon0;

    double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
    double xxpii, xypii, xzpii, yxpii, yypii, yzpii, zxpii, zypii, zzpii;
    double phidxp, phidyp, phidzp;
    double phidxpii, phidypii, phidzpii;
    double phidxpx, phidypx, phidzpx;
    double phidxpy, phidypy, phidzpy;
    double phidxpz, phidypz, phidzpz;

    double ep, epdx, epdy, epdz;

    double term0;
    double termx, termx0, termx1, termx0dx, termx1dx;
    double termy, termy0, termy1, termy0dy, termy1dy;
    double termz, termz0, termz1, termz0dz, termz1dz;

    double termiikk, termjjkk;

    double miijj;
    // *************  For anisotropy calculation **************

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
                if (ix == 0)
                {
                    ixm = 0;
                }
                if (ix == ndmx)
                {
                    ixp = ndmx;
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
                    izp = 0;
                }
                if (iz == 0)
                {
                    izm = ndmz;
                }
                double dF = 0.0;
                for (n1 = 1; n1 <= phiNum[ix][iy][iz]; n1++)
                {
                    ii = phiIdx[n1][ix][iy][iz];

                    phidxii = (phi[ii][ixp][iy][iz] - phi[ii][ixm][iy][iz]) / 2.0 / dx;
                    phidyii = (phi[ii][ix][iyp][iz] - phi[ii][ix][iym][iz]) / 2.0 / dx;
                    phidzii = (phi[ii][ix][iy][izp] - phi[ii][ix][iy][izm]) / 2.0 / dx;
                    phiabsii = phidxii * phidxii + phidyii * phidyii + phidzii * phidzii;

                    pddtt = 0.0;
                    for (n2 = 1; n2 <= phiNum[ix][iy][iz]; n2++)
                    {
                        jj = phiIdx[n2][ix][iy][iz];
                        intsum = 0.0;
                        for (n3 = 1; n3 <= phiNum[ix][iy][iz]; n3++)
                        {
                            kk = phiIdx[n3][ix][iy][iz];

                            // *************  For anisotropy calculation **************
                            phidx = (phi[kk][ixp][iy][iz] - phi[kk][ixm][iy][iz]) / 2.0 / dx;
                            phidy = (phi[kk][ix][iyp][iz] - phi[kk][ix][iym][iz]) / 2.0 / dx;
                            phidz = (phi[kk][ix][iy][izp] - phi[kk][ix][iy][izm]) / 2.0 / dx;

                            phidxx = (phi[kk][ixp][iy][iz] + phi[kk][ixm][iy][iz] - 2.0 * phi[kk][ix][iy][iz]) / dx / dx;
                            phidyy = (phi[kk][ix][iyp][iz] + phi[kk][ix][iym][iz] - 2.0 * phi[kk][ix][iy][iz]) / dx / dx;
                            phidzz = (phi[kk][ix][iy][izp] + phi[kk][ix][iy][izm] - 2.0 * phi[kk][ix][iy][iz]) / dx / dx;

                            phidxy = (phi[kk][ixp][iyp][iz] + phi[kk][ixm][iym][iz] - phi[kk][ixm][iyp][iz] - phi[kk][ixp][iym][iz]) / 4.0 / dx / dx;
                            phidxz = (phi[kk][ixp][iy][izp] + phi[kk][ixm][iy][izm] - phi[kk][ixm][iy][izp] - phi[kk][ixp][iy][izm]) / 4.0 / dx / dx;
                            phidyz = (phi[kk][ix][iyp][izp] + phi[kk][ix][iym][izm] - phi[kk][ix][iym][izp] - phi[kk][ix][iyp][izm]) / 4.0 / dx / dx;

                            phiabs = phidx * phidx + phidy * phidy + phidz * phidz;

                            if (anij[ii][kk] != 0.0 && phiabs != 0.0)
                            {
                                epsilon0 = sqrt(aij[ii][kk]);

                                th = thij[ii][kk];
                                vp = vpij[ii][kk];
                                eta = etaij[ii][kk];

                                xxp = cos(th) * cos(vp);
                                yxp = sin(th) * cos(vp);
                                zxp = sin(vp);
                                xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                zyp = cos(vp) * sin(eta);
                                xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                zzp = cos(eta) * cos(vp);

                                phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                termiikk = term0 + termx + termy + termz;
                            }
                            else
                            {
                                termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                            }

                            if (anij[jj][kk] != 0.0 && phiabs != 0.0)
                            {
                                epsilon0 = sqrt(aij[jj][kk]);

                                th = thij[jj][kk];
                                vp = vpij[jj][kk];
                                eta = etaij[jj][kk];

                                xxp = cos(th) * cos(vp);
                                yxp = sin(th) * cos(vp);
                                zxp = sin(vp);
                                xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                zyp = cos(vp) * sin(eta);
                                xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                zzp = cos(eta) * cos(vp);

                                phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                termjjkk = term0 + termx + termy + termz;
                            }
                            else
                            {
                                termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                            }
                            // ************************************************************

                            intsum += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][ix][iy][iz];
                        } // kk

                        thii = thij[ii][jj];
                        vpii = vpij[ii][jj];
                        etaii = etaij[ii][jj];

                        xxpii = cos(thii) * cos(vpii);
                        yxpii = sin(thii) * cos(vpii);
                        zxpii = sin(vpii);
                        xypii = -sin(thii) * cos(etaii) - cos(thii) * sin(vpii) * sin(etaii);
                        yypii = cos(thii) * cos(etaii) - sin(thii) * sin(vpii) * sin(etaii);
                        zypii = cos(vpii) * sin(etaii);
                        xzpii = sin(etaii) * sin(thii) - cos(etaii) * cos(thii) * sin(vpii);
                        yzpii = -sin(etaii) * cos(thii) - cos(etaii) * sin(thii) * sin(vpii);
                        zzpii = cos(etaii) * cos(vpii);

                        phidxpii = phidxii * xxpii + phidyii * yxpii + phidzii * zxpii;
                        phidypii = phidxii * xypii + phidyii * yypii + phidzii * zypii;
                        phidzpii = phidxii * xzpii + phidyii * yzpii + phidzii * zzpii;

                        if (anij[ii][jj] == 1.0 && phiabsii != 0.0)
                        {
                            miijj = mij[ii][jj] * (1.0 - 3.0 * astrem + 4.0 * astrem * (pow(phidxpii, 4.0) + pow(phidypii, 4.0) + pow(phidzpii, 4.0)) / pow(phiabsii, 2.0));
                        }
                        else
                        {
                            miijj = mij[ii][jj];
                        }
                        if (ii == 1 && jj == 0)
                        {
                            dF = -(temp[ix][iy][iz] - Tm) * dH / Tm;
                        }
                        else if (ii == 0 && jj == 1)
                        {
                            dF = (temp[ix][iy][iz] - Tm) * dH / Tm;
                        }
                        pddtt += -2.0 * miijj / double(phiNum[ix][iy][iz]) * (intsum - 8.0 / PI * dF * sqrt(phi[ii][ix][iy][iz] * phi[jj][ix][iy][iz]));
                    } // jj
                    phi2[ii][ix][iy][iz] = phi[ii][ix][iy][iz] + pddtt * dtime;
                    if (phi2[ii][ix][iy][iz] >= 1.0)
                    {
                        phi2[ii][ix][iy][iz] = 1.0;
                    }
                    if (phi2[ii][ix][iy][iz] <= 0.0)
                    {
                        phi2[ii][ix][iy][iz] = 0.0;
                    }
                    // termperature increase from release of latent heat
                    if (ii == 1)
                    {
                        temp[ix][iy][iz] += pddtt * dtime * dH * 1.5 / sph_s;
                    }
                } // ii
            }     // j
        }         // i
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
                    psum += phi2[kk][ix][iy][iz];
                }
                for (kk = 0; kk <= nm; kk++)
                {
                    phi2[kk][ix][iy][iz] = phi2[kk][ix][iy][iz] / psum;
                }
            }
        }
    }

    for (ix = start; ix <= end; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                for (kk = 0; kk <= nm; kk++)
                {
                    phi[kk][ix][iy][iz] = phi2[kk][ix][iy][iz];
                }
            }
        }
    }
}

void ComputeTemperature(double ***temp, double ***temp2, double ****phi, double T_right, double T_left,
                        int start, int end, int ndmx, int ndmy, int ndmz, double dtime, double dx,
                        int ix, int iy, int iz, int ixp, int ixm, int iyp, int iym, int izp, int izm,
                        double Dts, double Dtl)
{
    double Tddtt = 0.0;
    for (ix = start; ix <= end; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                double tempip = 0.0;
                double tempim = 0.0;
                ixp = ix + 1;
                ixm = ix - 1;
                iyp = iy + 1;
                iym = iy - 1;
                izp = iz + 1;
                izm = iz - 1;
                if (ix == ndmx)
                {
                    tempip = T_right;
                }
                else
                {
                    tempip = temp[ixp][iy][iz];
                }
                if (ix == 0)
                {
                    tempim = T_left;
                }
                else
                {
                    tempim = temp[ixm][iy][iz];
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
                    izp = 0;
                }
                if (iz == 0)
                {
                    izm = ndmz;
                }

                Tddtt = (Dtl * phi[0][ix][iy][iz] + Dts * phi[1][ix][iy][iz]) * (tempip + tempim + temp[ix][iyp][iz] + temp[ix][iym][iz] + temp[ix][iy][izp] + temp[ix][iy][izm] - 6.0 * temp[ix][iy][iz]) / dx / dx;
                temp2[ix][iy][iz] = temp[ix][iy][iz] + Tddtt * dtime;
            }
        }
    }

    for (ix = start; ix <= end; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                temp[ix][iy][iz] = temp2[ix][iy][iz];
            }
        }
    }
}

void MovingFrame(double ****phi, double ***temp,
                 double &T_right, double &T_left, double Tg,
                 int &intpos, int &curst, int &frapass,
                 int ndmx, int ndmy, int ndmz, int NDX, int NDY, int NDZ, double dx,
                 int ix, int iy, int iz, int ii, int nm, int istep,
                 double &int_temp)
{
    double sumplane = 0.0;
    int hasS = 0;
    int allS = 0;
    // check if the bottom is solid
    for (ix = 0; ix <= ndmx; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            if (phi[0][0][iy][iz] != 0.0)
            {
                hasS = 0;
            }
            sumplane += phi[0][0][iy][iz];
            if ((sumplane == 0.0) && (iy == ndmy) && (iz == ndmz))
            {
                hasS = 1;
            }
        }
    }
    // search interface front
    if (hasS == 1)
    {
        allS = 1;
        for (ix = 0; ix <= ndmx; ix++)
        {
            if (allS == 0)
            {
                intpos = ix - 1;
                int_temp = temp[ix][0][0];
                break;
            }
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    if (phi[0][ix][iy][iz] > 0.0)
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
    }

    if (intpos > NDX / 4)
    {
        frapass += 1;
        curst = istep;
        for (ix = 0; ix <= ndmx - 1; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    for (ii = 0; ii <= nm; ii++)
                    {
                        phi[ii][ix][iy][iz] = phi[ii][ix + 1][iy][iz];
                    }
                    temp[ix][iy][iz] = temp[ix + 1][iy][iz];
                }
            }
        }
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                for (ii = 0; ii <= nm; ii++)
                {
                    phi[ii][ndmx][iy][iz] = phi[ii][ndmx - 1][iy][iz];
                }
                temp[ndmx][iy][iz] = temp[ndmx - 1][iy][iz] + Tg * dx;
            }
        }
        T_left += Tg * dx;
        T_right += Tg * dx;
    }
}

void SavaData1D(double ****phi, double ***temp, int istep,
                int NDX, int NDY, int NDZ,
                int ndmx, int ndmy, int ndmz,
                int nm, int i, int j, int k)
{
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "data/phi/1d%d.csv", istep);
    stream = fopen(buffer, "a");

    for (i = 0; i <= ndmx; i++)
    {
        fprintf(stream, "%e\n", phi[1][i][0][0]);
    }
    fclose(stream);

    FILE *streamt;
    char buffert[30];
    sprintf(buffert, "data/temp/1d%d.csv", istep);
    streamt = fopen(buffert, "a");

    for (i = 0; i <= ndmx; i++)
    {
        fprintf(streamt, "%e\n", temp[i][0][0]);
    }
    fclose(streamt);
}