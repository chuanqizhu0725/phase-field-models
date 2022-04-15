
#include "header.h"

double ****phi, ****phi2;
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

    PhaseParameters(delta, gamma0, mobi, temp,
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
                           anij, thij, vpij, etaij, astre,
                           start, end, ndmx, ndmy, ndmz, dtime,
                           ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                           ii, jj, kk, n1, n2, n3, nm,
                           pddtt, intsum, psum);

#pragma omp barrier
        istep++;
        if (istep < nstep)
        {
            goto start;
        }
    end:;
    }

    SavaData3D(phi,
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

    return 0;
}
