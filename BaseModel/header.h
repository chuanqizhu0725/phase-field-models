#include <stdio.h>
#include <stdlib.h>
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
                     double A0, double W0, double M0, double S0,
                     double **aij, double **wij, double **mij, double **sij,
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
            sij[ni][nj] = 0.0;
            anij[ni][nj] = 0.0;
            thij[ni][nj] = 0.0;
            vpij[ni][nj] = 0.0;
            etaij[ni][nj] = 0.0;
            if ((ni == 0) || (nj == 0))
            {
                sij[ni][nj] = S0;
            }
            if (ni < nj)
            {
                sij[ni][nj] = -sij[ni][nj];
            }
            if (ni == nj)
            {
                wij[ni][nj] = 0.0;
                aij[ni][nj] = 0.0;
                mij[ni][nj] = 0.0;
                sij[ni][nj] = 0.0;
            }
        }
    }
    anij[1][0] = 1.0;
    anij[0][1] = 1.0;
}

void RandomSeeds(double ****phi,
                 int NDX, int NDY, int NDZ, int ndmx, int ndmy, int ndmz,
                 int i, int j, int k, int ni, int nj, int nm,
                 double r0, double r, int xx0, int yy0, int zz0)
{
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                phi[0][i][j][k] = 1.0;
                for (ni = 1; ni <= nm; ni++)
                {
                    phi[ni][i][j][k] = 0.0;
                }
            }
        }
    }

    r0 = 10.0;
    for (ni = 1; ni <= nm; ni++)
    {
        xx0 = rand() % NDX;
        yy0 = rand() % NDY;
        zz0 = rand() % NDZ;
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    r = sqrt((i - xx0) * (i - xx0) + (j - yy0) * (j - yy0) + (k - zz0) * (k - zz0));
                    if (r <= r0)
                    {
                        phi[ni][i][j][k] = 1.0;
                        phi[0][i][j][k] = 0.0;
                        for (nj = 1; nj <= nm; nj++)
                        {
                            if (nj != ni)
                            {
                                phi[nj][i][j][k] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
}

void TemperatureGradient(double ***temp, double temp0, double ****phi,
                         int ndmx, int ndmy, int ndmz,
                         int i, int j, int k)
{
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                if (phi[0][i][j][k] == 0.0)
                {
                    temp[i][j][k] = temp0;
                }
                else
                {
                    temp[i][j][k] = temp0 * 0.8;
                }
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
                r = sqrt((i - xx0) * (i - xx0) + (j - yy0) * (j - yy0) + (k - zz0) * (k - zz0));
                if (r < r0)
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

void CenterSeedBinary(double ****phi, double ***cont, double ****conp, double c0, double ***temp,
                      int NDX, int NDY, int NDZ, int ndmx, int ndmy, int ndmz,
                      int i, int j, int k, int nm,
                      double cl, double r0, double r, int xx0, int yy0, int zz0,
                      double Te, double ce, double ml1, double kap1, double ml2, double kap2)
{
    double sum = 0.0;
    double c1e = 0.0;
    double c01e = 0.0;
    double c2e = 0.0;
    double c02e = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                double c1e = (ce + (temp[i][j][k] - Te) / ml1) * kap1;
                double c01e = ce + (temp[i][j][k] - Te) / ml1;
                double c2e = 1.0 - (1.0 - (ce + (temp[i][j][k] - Te) / ml2)) * kap2;
                double c02e = ce + (temp[i][j][k] - Te) / ml2;
                r = sqrt((i - xx0) * (i - xx0) + (j - yy0) * (j - yy0) + (k - zz0) * (k - zz0));
                if (r < r0)
                {
                    phi[0][i][j][k] = 0.0;
                    conp[0][i][j][k] = c02e;
                    phi[nm][i][j][k] = 1.0;
                    conp[nm][i][j][k] = c2e;
                }
                else
                {
                    phi[0][i][j][k] = 1.0;
                    conp[0][i][j][k] = cl;
                    phi[nm][i][j][k] = 0.0;
                    conp[nm][i][j][k] = c2e;
                }
                cont[i][j][k] = conp[nm][i][j][k] * phi[nm][i][j][k] + conp[0][i][j][k] * phi[0][i][j][k];
                sum += cont[i][j][k];
            }
        }
    }
    c0 = sum / NDX / NDY / NDZ;
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
                        double ****conp, double ***temp, double Te, double ce, double ml1, double ml2)
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
                            dF = ((conp[0][ix][iy][iz] - ce) * ml1 + Te - temp[ix][iy][iz]) * sij[1][0];
                        }
                        else if (ii == 0 && jj == 1)
                        {
                            dF = ((conp[0][ix][iy][iz] - ce) * ml1 + Te - temp[ix][iy][iz]) * sij[0][1];
                        }
                        else if (ii == 2 && jj == 0)
                        {
                            dF = ((conp[0][ix][iy][iz] - ce) * ml2 + Te - temp[ix][iy][iz]) * sij[2][0];
                        }
                        else if (ii == 0 && jj == 2)
                        {
                            dF = -((conp[0][ix][iy][iz] - ce) * ml2 + Te - temp[ix][iy][iz]) * sij[2][0];
                        }
                        else
                        {
                            dF = 0.0;
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

void SolutePartition(double ****conp, double ***cont, double ***cont2, double ****phi, double ***temp,
                     int start, int end, int ndmx, int ndmy, int ndmz, double dtime, double dx,
                     int ix, int iy, int iz, int ixp, int ixm, int iyp, int iym, int izp, int izm,
                     double Te, double ce, double ml1, double kap1, double ml2, double kap2)
{
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
                    izp = 0;
                }
                if (iz == 0)
                {
                    izm = ndmz;
                }

                // Calculate equilirium concentration
                double c1e = (ce + (temp[ix][iy][iz] - Te) / ml1) * kap1;
                double c01e = ce + (temp[ix][iy][iz] - Te) / ml1;
                double c2e = 1.0 - (1.0 - (ce + (temp[ix][iy][iz] - Te) / ml2)) * kap2;
                double c02e = ce + (temp[ix][iy][iz] - Te) / ml2;

                if (phi[0][ix][iy][iz] == 0.0)
                {
                    if (phi[1][ix][iy][iz] == 1.0)
                    {
                        conp[1][ix][iy][iz] = cont[ix][iy][iz];
                    }
                    else
                    {
                        conp[1][ix][iy][iz] = c1e;
                    }
                    if (phi[2][ix][iy][iz] == 1.0)
                    {
                        conp[2][ix][iy][iz] = cont[ix][iy][iz];
                    }
                    else
                    {
                        conp[2][ix][iy][iz] = c2e;
                    }
                    // Correct abnormal calculation at solid edge
                    if ((phi[0][ixp][iy][iz] > 0.0) || (phi[0][ixm][iy][iz] > 0.0) || (phi[0][ix][iyp][iz] > 0.0) || (phi[0][ix][iym][iz] > 0.0) || (phi[0][ix][iy][izp] > 0.0) || (phi[0][ix][iy][izm] > 0.0))
                    {
                        conp[1][ix][iy][iz] = c1e;
                        conp[2][ix][iy][iz] = c2e;
                    }
                    conp[0][ix][iy][iz] = c01e * phi[1][ix][iy][iz] + c02e * phi[2][ix][iy][iz];
                }
                else if (phi[0][ix][iy][iz] > 0.0 && phi[0][ix][iy][iz] < 1.0)
                {
                    conp[1][ix][iy][iz] = c1e;
                    conp[2][ix][iy][iz] = c2e;
                    conp[0][ix][iy][iz] = (cont[ix][iy][iz] - conp[1][ix][iy][iz] * phi[1][ix][iy][iz] - conp[2][ix][iy][iz] * phi[2][ix][iy][iz]) / phi[0][ix][iy][iz];
                    // Correct abnormal calculation at liquid edge
                    if (phi[0][ix][iy][iz] < 0.05)
                    {
                        conp[0][ix][iy][iz] = (c01e * phi[1][ix][iy][iz] + c02e * phi[2][ix][iy][iz]) / (phi[1][ix][iy][iz] + phi[2][ix][iy][iz]);
                    }
                    if (conp[0][ix][iy][iz] > 1.0)
                    {
                        conp[0][ix][iy][iz] = 1.0;
                    }
                    if (conp[0][ix][iy][iz] < 0.0)
                    {
                        conp[0][ix][iy][iz] = 0.0;
                    }
                }
                else if (phi[0][ix][iy][iz] == 1.0)
                {
                    conp[1][ix][iy][iz] = c1e;
                    conp[2][ix][iy][iz] = c2e;
                    conp[0][ix][iy][iz] = cont[ix][iy][iz];
                }
                cont[ix][iy][iz] = conp[1][ix][iy][iz] * phi[1][ix][iy][iz] + conp[2][ix][iy][iz] * phi[2][ix][iy][iz] + conp[0][ix][iy][iz] * phi[0][ix][iy][iz];
                if (cont[ix][iy][iz] > 1.0)
                {
                    cont[ix][iy][iz] = 1.0;
                }
                if (cont[ix][iy][iz] < 0.0)
                {
                    cont[ix][iy][iz] = 0.0;
                }
            }
        }
    }
}

void ComputeConcentration(double ****conp, double ***cont, double ***cont2, double ****phi,
                          int start, int end, int ndmx, int ndmy, int ndmz, double dtime, double dx,
                          int ix, int iy, int iz, int ixp, int ixm, int iyp, int iym, int izp, int izm,
                          int ii, int nm, double cddtt,
                          double Dcs, double Dcl)
{
    double sumcs = 0.0;
    double sumcl = 0.0;
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
                    izp = 0;
                }
                if (iz == 0)
                {
                    izm = ndmz;
                }
                //拡散方程式内における微分計算
                for (ii = 1; ii < nm; ii++)
                {
                    sumcs = 0.25 * ((phi[ii][ixp][iy][iz] - phi[ii][ixm][iy][iz]) * (conp[ii][ixp][iy][iz] - conp[ii][ixm][iy][iz]) + (phi[ii][ix][iyp][iz] - phi[ii][ix][iym][iz]) * (conp[ii][ix][iyp][iz] - conp[ii][ix][iym][iz]) + (phi[ii][ix][iy][izp] - phi[ii][ix][iy][izm]) * (conp[ii][ix][iy][izp] - conp[ii][ix][iy][izm])) / dx / dx +
                            phi[ii][ix][iy][iz] * (conp[ii][ixp][iy][iz] + conp[ii][ixm][iy][iz] + conp[ii][ix][iyp][iz] + conp[ii][ix][iym][iz] + conp[ii][ix][iy][izp] + conp[ii][ix][iy][izm] - 6.0 * conp[ii][ix][iy][iz]) / dx / dx;
                }
                sumcl = 0.25 * ((phi[0][ixp][iy][iz] - phi[0][ixm][iy][iz]) * (conp[0][ixp][iy][iz] - conp[0][ixm][iy][iz]) + (phi[0][ix][iyp][iz] - phi[0][ix][iym][iz]) * (conp[0][ix][iyp][iz] - conp[0][ix][iym][iz]) + (phi[0][ix][iy][izp] - phi[0][ix][iy][izm]) * (conp[0][ix][iy][izp] - conp[0][ix][iy][izm])) / dx / dx +
                        phi[0][ix][iy][iz] * (conp[0][ixp][iy][iz] + conp[0][ixm][iy][iz] + conp[0][ix][iyp][iz] + conp[0][ix][iym][iz] + conp[0][ix][iy][izp] + conp[0][ix][iy][izm] - 6.0 * conp[0][ix][iy][iz]) / dx / dx;
                cddtt = Dcs * sumcs + Dcl * sumcl;
                cont2[ix][iy][iz] = cont[ix][iy][iz] + cddtt * dtime;
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
                cont[ix][iy][iz] = cont2[ix][iy][iz];
            }
        }
    }
}

void ComputeTemperature(double ***temp, double ***temp2, double ****phi,
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
                    izp = 0;
                }
                if (iz == 0)
                {
                    izm = ndmz;
                }
                Tddtt = (Dtl * phi[0][ix][iy][iz] + Dts * phi[1][ix][iy][iz]) * (temp[ixp][iy][iz] + temp[ixm][iy][iz] + temp[ix][iyp][iz] + temp[ix][iym][iz] + temp[ix][iy][izp] + temp[ix][iy][izm] - 6.0 * temp[ix][iy][iz]) / dx / dx;
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
    for (iy = 0; iy <= ndmy; iy++)
    {
        for (iz = 0; iz <= ndmz; iz++)
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
                int_temp = temp[intpos][0][0];
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

    if (intpos > (NDX / 4))
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

void SavaData3D(double ****phi, int istep,
                int NDX, int NDY, int NDZ,
                int ndmx, int ndmy, int ndmz,
                int nm, int i, int j, int k)
{
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "data/phi/3d%d.vtk", istep);
    stream = fopen(buffer, "a");

    fprintf(stream, "# vtk DataFile Version 1.0\n");
    fprintf(stream, "phi.vtk\n");
    fprintf(stream, "ASCII\n");
    fprintf(stream, "DATASET STRUCTURED_POINTS\n");
    fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
    fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
    fprintf(stream, "\n");
    fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
    fprintf(stream, "SCALARS scalars float\n");
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (k = 0; k <= ndmz; k++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (i = 0; i <= ndmx; i++)
            {
                fprintf(stream, "%e\n", phi[1][i][j][k]);
            }
        }
    }
    fclose(stream);
}

void SavaData2D(double ****phi, int istep,
                int NDX, int NDY, int NDZ,
                int ndmx, int ndmy, int ndmz,
                int nm, int i, int j, int k)
{
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "2d%d.csv", istep);
    stream = fopen(buffer, "a");

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            fprintf(stream, "%e\n", phi[1][i][j][0]);
        }
    }
    fclose(stream);
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
        fprintf(stream, "%e\n", phi[2][i][0][0]);
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

void SavaData1DBinary(double ****phi, double ***cont, int istep,
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
        fprintf(stream, "%e\n", phi[2][i][0][0]);
    }
    fclose(stream);

    FILE *streamc;
    char bufferc[30];
    sprintf(bufferc, "data/con/1d%d.csv", istep);
    streamc = fopen(bufferc, "a");

    for (i = 0; i <= ndmx; i++)
    {
        fprintf(streamc, "%e\n", cont[i][0][0]);
    }
    fclose(streamc);
}