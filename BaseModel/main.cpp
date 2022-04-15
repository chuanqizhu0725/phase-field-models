
#include "header.h"

int N = 2;
int NTH = 8;
int NDX = 64;
int NDY = 64;
int NDZ = 64;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;

int nstep = 401;
int pstep = 100;
double dx0 = 2.0e-8;
double dx = 1.0;
double dtime = 5.0;
double temp = 1000.0;
double vm0 = 7.0e-6;
double delta = 7.0;
double mobi = 1.0;
double astre = 0.05;
double astrem = 0.5;
double gamma0 = 0.5 * vm0 / RR / temp / dx0;

double Te = 0.0;
double ce = 0.5;
double ml1 = -10.0;
double kap1 = 0.2;
double ml2 = 10.0;
double kap2 = 0.2;

int i, j, k, ni, nj;
int xx0, yy0, zz0;
double r0, r;
double M0, W0, A0, F0;

CImg<unsigned char> phi_fldxz(NDX, NDZ, 1, 3), phi_fldxy(NDX, NDY, 1, 3);
char outFilePhi_xz[64], outFilePhi_xy[64];

double ****phi, ****phi2;
double ****conp, ***cont, ***cont2;
int ***phiNum, ****phiIdx;
double **aij, **wij, **mij, **fij;
double **anij, **thij, **vpij, **etaij;

int main(int argc, char *argv[])
{
    phi = new double ***[N];
    phi2 = new double ***[N];
    for (ni = 0; ni <= nm; ni++)
    {
        phi[ni] = new double **[NDX];
        phi2[ni] = new double **[NDX];
        for (i = 0; i <= ndmx; i++)
        {
            phi[ni][i] = new double *[NDY];
            phi2[ni][i] = new double *[NDY];
            for (j = 0; j <= ndmy; j++)
            {
                phi[ni][i][j] = new double[NDZ];
                phi2[ni][i][j] = new double[NDZ];
            }
        }
    }

    phiIdx = new int ***[N + 1];
    for (ni = 0; ni <= N; ni++)
    {
        phiIdx[ni] = new int **[NDX];
        for (i = 0; i <= ndmx; i++)
        {
            phiIdx[ni][i] = new int *[NDY];
            for (j = 0; j <= ndmy; j++)
            {
                phiIdx[ni][i][j] = new int[NDZ];
            }
        }
    }

    phiNum = new int **[NDX];
    for (i = 0; i <= ndmx; i++)
    {
        phiNum[i] = new int *[NDY];

        for (j = 0; j <= ndmy; j++)
        {
            phiNum[i][j] = new int[NDZ];
        }
    }

    aij = new double *[N];
    wij = new double *[N];
    mij = new double *[N];
    fij = new double *[N];
    anij = new double *[N];
    thij = new double *[N];
    vpij = new double *[N];
    etaij = new double *[N];
    for (ni = 0; ni <= nm; ni++)
    {
        aij[ni] = new double[N];
        wij[ni] = new double[N];
        mij[ni] = new double[N];
        fij[ni] = new double[N];
        anij[ni] = new double[N];
        thij[ni] = new double[N];
        vpij[ni] = new double[N];
        etaij[ni] = new double[N];
    }

    PhaseProperties(delta, gamma0, mobi, temp,
                    A0, W0, M0, F0,
                    aij, wij, mij, fij,
                    anij, thij, vpij, etaij,
                    ni, nj, nm);

    // RandomSeeds(phi,
    //             NDX, NDY, NDZ, ndmx, ndmy, ndmz,
    //             i, j, k, ni, nj, nm,
    //             r0, r, xx0, yy0, zz0);

    CenterSeed(phi,
               NDX, NDY, NDZ, ndmx, ndmy, ndmz,
               i, j, k, nm,
               NDX / 8, r, NDX / 2, NDY / 2, NDZ / 2);

    int rows = NDX / NTH;

#pragma omp parallel num_threads(NTH)
    {
        int th_id, offset, start, end, istep;
        int ix, iy, iz, ixp, ixm, iyp, iym, izp, izm;
        int ii, jj, kk, phinum;
        int n1, n2, n3;
        double intsum, pddtt, psum;
        double cddtt, sumcs, sumcl;

        th_id = omp_get_thread_num();
        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;
        istep = 0;

    start:;

        CollectPhaseFields(phi, phiNum, phiIdx,
                           start, end, ndmx, ndmy, ndmz,
                           ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                           phinum, nm);

        ComputePhaseFields(phi, phi2, phiNum, phiIdx,
                           aij, wij, mij, fij,
                           anij, thij, vpij, etaij, astre, astrem,
                           start, end, ndmx, ndmy, ndmz, dtime, dx,
                           ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                           ii, jj, kk, n1, n2, n3, nm,
                           pddtt, intsum, psum);

#pragma omp barrier

        CollectPhaseFields(phi, phiNum, phiIdx,
                           start, end, ndmx, ndmy, ndmz,
                           ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                           phinum, nm);

        SolutePartition(conp, cont, cont2, phi, temp,
                        start, end, ndmx, ndmy, ndmz, dtime, dx,
                        ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                        Te, ce, ml1, kap1, ml2, kap2);

#pragma omp barrier

        ComputeConcentration();

#pragma omp barrier
        if ((istep % pstep == 0) && (th_id == 0))
        {
            // ****** XZ *******
            cimg_forXY(phi_fldxz, x, z)
            {
                phi_fldxz(x, z, 0) = 255. * (phi[1][x][NDY / 2][z]); // red
                phi_fldxz(x, z, 1) = 255. * (phi[1][x][NDY / 2][z]); // green
                phi_fldxz(x, z, 2) = 255. * (phi[1][x][NDY / 2][z]); // blue
            }
            sprintf(outFilePhi_xz, "figures/phi/2dxz%d.png", istep);
            phi_fldxz.save_jpeg(outFilePhi_xz);

            SavaData3D(phi, istep,
                       NDX, NDY, NDZ,
                       ndmx, ndmy, ndmz,
                       nm, i, j, k);

            // SavaData2D(phi,
            //            NDX, NDY, NDZ,
            //            ndmx, ndmy, ndmz,
            //            nm, i, j, k);

            // SavaData1D(phi,
            //            NDX, NDY, NDZ,
            //            ndmx, ndmy, ndmz,
            //            nm, i, j, k);
        }
        istep++;
        if (istep < nstep)
        {
            goto start;
        }
    end:;
    }
    return 0;
}
