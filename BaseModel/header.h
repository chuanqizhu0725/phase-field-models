#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <omp.h>

using namespace std;

#define PI 3.141592
#define RR 8.3145

int N = 3;
int NTH = 8;
int NDX = 128;
int NDY = 1;
int NDZ = 1;
int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;

int nstep = 1001;
double dx = 2.0e-8;
double dtime = 5.0;
double temp = 1000.0;
double vm0 = 7.0e-6;
double delta = 7.0;
double mobi = 1.0;
double gamma0 = 0.5 * vm0 / RR / temp / dx;

int i, j, k, ni, nj;
int xx0, yy0, zz0;
double r0, r;
double M0, W0, A0, F0;

double ****phi, ****phi2;
int ***phiNum, ****phiIdx;
double **aij, **wij, **mij, **fij;

void PhaseParameters(double delta, double gamma0, double mobi, double temp,
                     double A0, double W0, double M0, double F0,
                     double **aij, double **wij, double **mij, double **fij,
                     int ni, int nj, int nm)
{
    A0 = 8.0 * delta * gamma0 / PI / PI;
    W0 = 4.0 * gamma0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);
    F0 = 80.0 / RR / temp;

    for (ni = 0; ni <= nm; ni++)
    {
        for (nj = 0; nj <= nm; nj++)
        {
            wij[ni][nj] = W0;
            aij[ni][nj] = A0;
            mij[ni][nj] = M0;
            fij[ni][nj] = 0.0;
            if ((ni == 0) || (nj == 0))
            {
                fij[ni][nj] = F0;
            }
            if (ni < nj)
            {
                fij[ni][nj] = -fij[ni][nj];
            }
            if (ni == nj)
            {
                wij[ni][nj] = 0.0;
                aij[ni][nj] = 0.0;
                mij[ni][nj] = 0.0;
                fij[ni][nj] = 0.0;
            }
        }
    }
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
                        double **aij, double **wij, double **mij, double **fij,
                        int start, int end, int ndmx, int ndmy, int ndmz, double dtime,
                        int ix, int iy, int iz, int ixp, int ixm, int iyp, int iym, int izp, int izm,
                        int ii, int jj, int kk, int n1, int n2, int n3, int nm,
                        double pddtt, double intsum, double psum)
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

                for (n1 = 1; n1 <= phiNum[ix][iy][iz]; n1++)
                {
                    ii = phiIdx[n1][ix][iy][iz];
                    pddtt = 0.0;
                    for (n2 = 1; n2 <= phiNum[ix][iy][iz]; n2++)
                    {
                        jj = phiIdx[n2][ix][iy][iz];
                        intsum = 0.0;
                        for (n3 = 1; n3 <= phiNum[ix][iy][iz]; n3++)
                        {
                            kk = phiIdx[n3][ix][iy][iz];
                            intsum += 0.5 * (aij[ii][kk] - aij[jj][kk]) * (phi[kk][ixp][iy][iz] + phi[kk][ixm][iy][iz] + phi[kk][ix][iyp][iz] + phi[kk][ix][iym][iz] + phi[kk][ix][iy][izp] + phi[kk][ix][iy][izm] - 6.0 * phi[kk][ix][iy][iz]) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][ix][iy][iz]; //[式(4.31)の一部]
                        }
                        pddtt += -2.0 * mij[ii][jj] / double(phiNum[ix][iy][iz]) * (intsum - 8.0 / PI * fij[ii][jj] * sqrt(phi[ii][ix][iy][iz] * phi[jj][ix][iy][iz]));
                    }
                    phi2[ii][ix][iy][iz] = phi[ii][ix][iy][iz] + pddtt * dtime;
                    if (phi2[ii][ix][iy][iz] >= 1.0)
                    {
                        phi2[ii][ix][iy][iz] = 1.0;
                    }
                    if (phi2[ii][ix][iy][iz] <= 0.0)
                    {
                        phi2[ii][ix][iy][iz] = 0.0;
                    }
                }
            } // j
        }     // i
    }

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

    for (kk = 0; kk <= nm; kk++)
    {
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    phi[kk][ix][iy][iz] = phi2[kk][ix][iy][iz];
                }
            }
        }
    }
}

void SavaData3D(double ****phi,
                int NDX, int NDY, int NDZ,
                int ndmx, int ndmy, int ndmz,
                int nm, int i, int j, int k)
{
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "3d.vtk");
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

void SavaData2D(double ****phi,
                int NDX, int NDY, int NDZ,
                int ndmx, int ndmy, int ndmz,
                int nm, int i, int j, int k)
{
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "2d.csv");
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

void SavaData1D(double ****phi,
                int NDX, int NDY, int NDZ,
                int ndmx, int ndmy, int ndmz,
                int nm, int i, int j, int k)
{
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "2d.csv");
    stream = fopen(buffer, "a");

    for (i = 0; i <= ndmx; i++)
    {
        fprintf(stream, "%e\n", phi[1][i][0][0]);
    }
    fclose(stream);
}