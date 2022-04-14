
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

#define N 2
#define ND 100
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;
int mid = ND / 4;

int nstep = 5000001;
int pstep = 100000;

double dx = 1.0e-5;
double dtime = 1.0e-6;
double gamma0 = 0.5;
double mobi = 1.0e-11;
double delta = 5.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);

double Tm = 1687.0;
double sph_s = 2.29e6;
double kap_s = 22;
double sph_l = 2.53e6;
double kap_l = 54;

double Cdt_s = kap_s / sph_s;
double Cdt_l = kap_l / sph_l;
double dH = 4.122e9;

double Tg = 8.0e3;
double Tv = 2.0e-4;
double Tr = Tg * Tv;

double T_left = 1674.0 - ND / 4 * dx * Tg;
double T_right = T_left + Tg * ND * dx;

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

double phi[N][ND], phi2[N][ND];
double temp[ND], temp2[ND];
double tempip, tempim;

int phinum;
int phiNum[ND];
int phiIdx[N + 1][ND];

int i, j, im, ip, k;
int ii, jj, kk;
int n1, n2, n3;
int istep;
int intpos, dist, curpos, prepos, frapass;
int curt, pret;
double int_vel;
double int_temp;

double F0, pddtt, sum1;
double termiikk, termjjkk;

double Tddtt;

void datasave(int step);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "temperature field stablity number is: " << dtime * Cdt_s / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * M0 * A0 << endl;

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

    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        if (i <= ND / 4)
        {
            phi[1][i] = 1.0;
            phi[0][i] = 0.0;
        }
        else
        {
            phi[1][i] = 0.0;
            phi[0][i] = 1.0;
        }
        temp[i] = T_left + i * dx * Tg;
    }

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        curpos = frapass + intpos;
        int_vel = double(curpos - prepos) * dx / double(curt - pret) / dtime;
        prepos = curpos;
        pret = curt;
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
        if (int_vel > 0.0)
        {
            cout << "interface velocity is: " << int_vel << endl;
        }
        cout << "interface temperature is " << temp[intpos] << endl;
    }

    for (i = 0; i <= ndm; i++)
    {
        ip = i + 1;
        im = i - 1;
        if (i == ndm)
        {
            ip = ndm - 1;
        }
        if (i == 0)
        {
            im = 1;
        }

        phinum = 0;
        for (ii = 0; ii <= nm; ii++)
        {
            if ((phi[ii][i] > 0.0) ||
                ((phi[ii][i] == 0.0) && (phi[ii][ip] > 0.0) || (phi[ii][im] > 0.0)))
            {
                phinum++;
                phiIdx[phinum][i] = ii;
            }
        }
        phiNum[i] = phinum;
    }

    // Evolution Equation of Phase Fields
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

        for (n1 = 1; n1 <= phiNum[i]; n1++)
        {
            ii = phiIdx[n1][i];
            pddtt = 0.0;
            for (n2 = 1; n2 <= phiNum[i]; n2++)
            {
                jj = phiIdx[n2][i];
                sum1 = 0.0;
                for (n3 = 1; n3 <= phiNum[i]; n3++)
                {
                    kk = phiIdx[n3][i];

                    termiikk = aij[ii][kk] * (phi[kk][ip] + phi[kk][im] - 2.0 * phi[kk][i]) / (dx * dx);

                    termjjkk = aij[jj][kk] * (phi[kk][ip] + phi[kk][im] - 2.0 * phi[kk][i]) / (dx * dx);

                    sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i];
                }
                if (ii == 1 && jj == 0)
                {
                    F0 = -(temp[i] - Tm) * dH / Tm;
                }
                else if (ii == 0 && jj == 1)
                {
                    F0 = (temp[i] - Tm) * dH / Tm;
                }

                pddtt += -2.0 * mij[ii][jj] / double(phiNum[i]) * (sum1 - 8.0 / PI * F0 * sqrt(phi[ii][i] * phi[jj][i]));
            }
            phi2[ii][i] = phi[ii][i] + pddtt * dtime;
            if (phi2[ii][i] >= 1.0)
            {
                phi2[ii][i] = 1.0;
            }
            if (phi2[ii][i] <= 0.0)
            {
                phi2[ii][i] = 0.0;
            }
            // termperature increase from release of latent heat
            if (ii == 1)
            {
                temp[i] += pddtt * dtime * dH / sph_s;
            }
        }
    } // i

    //
    for (i = 0; i <= ndm; i++)
    {
        sum1 = 0.0;
        for (k = 0; k <= nm; k++)
        {
            sum1 += phi2[k][i];
        }
        for (k = 0; k <= nm; k++)
        {
            phi2[k][i] = phi2[k][i] / sum1;
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (k = 0; k <= nm; k++)
        {
            phi[k][i] = phi2[k][i];
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        ip = i + 1;
        im = i - 1;
        // right boundary
        if (i == ndm)
        {
            tempip = T_right;
        }
        else
        {
            tempip = temp[ip];
        }
        // left boundary
        if (i == 0)
        {
            tempim = T_left;
        }
        else
        {
            tempim = temp[im];
        }
        Tddtt = (Cdt_l * phi[0][i] + Cdt_s * phi[1][i]) * (tempip + tempim - 2.0 * temp[i]) / dx / dx;
        temp2[i] = temp[i] + Tddtt * dtime;
    }

    for (i = 0; i <= ndm; i++)
    {
        temp[i] = temp2[i];
    }

    T_left -= Tr * dtime;
    T_right -= Tr * dtime;

    intpos = 0;
    for (i = 0; i <= ndm; i++)
    {
        if (phi[1][i] < 1.0)
        {
            intpos = i;
            int_temp = temp[i];
            break;
        }
    }

    if (intpos > mid)
    {
        frapass += 1;
        curt = istep;
        for (i = 0; i <= ndm - 1; i++)
        {
            for (k = 0; k <= nm; k++)
            {
                phi[k][i] = phi[k][i + 1];
            }
            temp[i] = temp[i + 1];
        }
        for (k = 0; k <= nm; k++)
        {
            phi[k][ndm] = phi[k][ndm - 1];
        }
        temp[ndm] = temp[ndm - 1] + Tg * dx;
        T_left += Tg * dx;
        T_right += Tg * dx;
    }

    istep = istep + 1;
    if (istep < nstep)
    {
        goto start;
    }

end:;
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

    FILE *streamt; //ストリームのポインタ設定
    char buffert[30];
    sprintf(buffert, "data/temp/1d%d.csv", step);
    streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        fprintf(streamt, "%e   ", temp[i]);
        fprintf(streamt, "\n");
    }
    fclose(streamt); //ファイルをクローズ
}
