
#include "header.h"

int N = 2;
int NTH = 1;
int NDX = 100;
int NDY = 1;
int NDZ = 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;

int nstep = 1000001;
int pstep = 100000;
double dx = 1.0e-5;
double dtime = 1.0e-6;
double delta = 5.0 * dx;
double mobi = 1.0e-11;
double astre = 0.00;
double astrem = 0.0;
double gamma0 = 0.5;
double Tm = 1687.0;
double sph_s = 2.29e6;
double kap_s = 22.0;
double sph_l = 2.53e6;
double kap_l = 54.0;
double Dts = kap_s / sph_s;
double Dtl = kap_l / sph_l;
double dH = 4.122e9;

double Tg = 8.0e3;
double Tv = .5e-4;
double Tr = Tg * Tv;

double temp0 = 1682.09;
double T_left = temp0 - NDX / 4 * dx * Tg;
double T_right = T_left + Tg * NDX * dx;

int i, j, k, ni, nj;
int xx0, yy0, zz0;
double r0, r, c0;
double M0, W0, A0;

int intpos, dist, curpos, prepos, frapass;
int curst, prest;
double int_vel;
double int_temp;

CImg<unsigned char> phi_fldxz(NDX, NDZ, 1, 3), phi_fldxy(NDX, NDY, 1, 3);
char outFilePhi_xz[64], outFilePhi_xy[64];

double ****phi, ****phi2;
double ****conp, ***cont, ***cont2, ***temp, ***temp2;
int ***phiNum, ****phiIdx;
double **aij, **wij, **mij, **sij;
double **anij, **thij, **vpij, **etaij;

int main(int argc, char *argv[])
{
    phi = new double ***[N];
    phi2 = new double ***[N];
    conp = new double ***[N];
    for (ni = 0; ni <= nm; ni++)
    {
        phi[ni] = new double **[NDX];
        phi2[ni] = new double **[NDX];
        conp[ni] = new double **[NDX];
        for (i = 0; i <= ndmx; i++)
        {
            phi[ni][i] = new double *[NDY];
            phi2[ni][i] = new double *[NDY];
            conp[ni][i] = new double *[NDY];
            for (j = 0; j <= ndmy; j++)
            {
                phi[ni][i][j] = new double[NDZ];
                phi2[ni][i][j] = new double[NDZ];
                conp[ni][i][j] = new double[NDZ];
            }
        }
    }

    cont = new double **[NDX];
    cont2 = new double **[NDX];
    temp = new double **[NDX];
    temp2 = new double **[NDX];
    for (i = 0; i <= ndmx; i++)
    {
        cont[i] = new double *[NDY];
        cont2[i] = new double *[NDY];
        temp[i] = new double *[NDY];
        temp2[i] = new double *[NDY];
        for (j = 0; j <= ndmy; j++)
        {
            cont[i][j] = new double[NDZ];
            cont2[i][j] = new double[NDZ];
            temp[i][j] = new double[NDZ];
            temp2[i][j] = new double[NDZ];
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
    sij = new double *[N];
    anij = new double *[N];
    thij = new double *[N];
    vpij = new double *[N];
    etaij = new double *[N];
    for (ni = 0; ni <= nm; ni++)
    {
        aij[ni] = new double[N];
        wij[ni] = new double[N];
        mij[ni] = new double[N];
        sij[ni] = new double[N];
        anij[ni] = new double[N];
        thij[ni] = new double[N];
        vpij[ni] = new double[N];
        etaij[ni] = new double[N];
    }

    PhaseProperties(delta, gamma0, mobi,
                    A0, W0, M0,
                    aij, wij, mij,
                    anij, thij, vpij, etaij,
                    ni, nj, nm);

    CenterSeed(phi,
               NDX, NDY, NDZ, ndmx, ndmy, ndmz,
               i, j, k, nm,
               NDX / 4, r, 0, NDY / 2, NDZ / 2);

    TemperatureGradient(temp, temp0, T_left, Tg,
                        ndmx, ndmy, ndmz, dx,
                        i, j, k);

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

#pragma omp barrier

        CollectPhaseFields(phi, phiNum, phiIdx,
                           start, end, ndmx, ndmy, ndmz,
                           ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                           phinum, nm);

        ComputePhaseFields(phi, phi2, phiNum, phiIdx,
                           aij, wij, mij, sij,
                           anij, thij, vpij, etaij,
                           astre, astrem,
                           start, end, ndmx, ndmy, ndmz, dtime, dx,
                           ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                           ii, jj, kk, n1, n2, n3, nm,
                           pddtt, intsum, psum,
                           conp, temp, dH, Tm, sph_s);

#pragma omp barrier

        ComputeTemperature(temp, temp2, phi, T_right, T_left,
                           start, end, ndmx, ndmy, ndmz, dtime, dx,
                           ix, iy, iz, ixp, ixm, iyp, iym, izp, izm,
                           Dts, Dtl);

#pragma omp barrier

        if (th_id == 0)
        {
            T_left -= Tr * dtime;
            T_right -= Tr * dtime;

            MovingFrame(phi, temp,
                        T_right, T_left, Tg,
                        intpos, curst, frapass,
                        ndmx, ndmy, ndmz, NDX, NDY, NDZ, dx,
                        i, j, k, ii, nm, istep,
                        int_temp);
        }

        if ((istep % pstep == 0) && (th_id == 0))
        {
            curpos = intpos + frapass;
            int_vel = double(curpos - prepos) * dx / double(curst - prest) / dtime;
            prepos = curpos;
            prest = curst;
            SavaData1D(phi, temp, istep,
                       NDX, NDY, NDZ,
                       ndmx, ndmy, ndmz,
                       nm, i, j, k);

            std::cout << "interface postion is: " << intpos << std::endl;
            std::cout << "interface temperature is: " << int_temp << std::endl;
            if (int_vel > 0.0 && int_vel < 1.0)
            {
                std::cout << "interface velocity is: " << int_vel << std::endl;
            }
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
