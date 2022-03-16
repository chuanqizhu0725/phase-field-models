// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/2/16
// Updated on 2022/03/07

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>

#define NDX 50 //差分計算における計算領域一辺の分割数
#define NDY 50 //差分計算における計算領域一辺の分割数
#define N 3    //考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)

int ndx = NDX;
int ndy = NDY;
int ndmx = NDX - 1;
int ndmy = NDY - 1;          //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int nm = N - 1, nmm = N - 2; //考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
double PI = 3.141592;        //π、計算カウント数
double RR = 8.3145;          //ガス定数

double aij[N][N]; //勾配エネルギー係数
double wij[N][N]; //ペナルティー項の係数
double mij[N][N]; //粒界の易動度
double fij[N][N]; //粒界移動の駆動力
int anij[N][N];
double thij[N][N];
int phinum;

int i, j, k, l, ii, jj, kk, ll, it; //整数
int ip, im, jp, jm;                 //整数
int n1, n2, n3;                     //整数

int istep = 0;
// int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
int nstep;               //計算カウント数の最大値（計算終了カウント）
double dtime, L, dx;     // L計算領域の一辺の長さ(nm), 差分プロック１辺の長さ(m)
double M0;               //粒界の易動度
double W0;               //ペナルティー項の係数
double A0;               //勾配エネルギー係数
double F0;               //粒界移動の駆動力
double temp;             //温度
double sum1, sum2, sum3; //各種の和の作業変数
double pddtt;            //フェーズフィールドの時間変化率

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;    //モル体積

double astre;
double theta, theta0;
double thetax, thetay;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidxx, phidyy, phidxy;
double phiabs;
double phidxp, phidyp;
double phidxpx, phidypx, phidxpy, phidypy;
double ep, epdx, epdy;
double term0;
double termx, termx0, termx1, termx0dx, termx1dx;
double termy, termy0, termy1, termy0dy, termy1dy;

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    nstep = 701;
    dtime = 5.0;
    temp = 1000.0;
    L = 2000.0;
    vm0 = 7.0e-6;
    delta = 7.0;
    mobi = 1.0;

    dx = L / (double)NDX * 1.0e-9;       //差分プロック１辺の長さ(m)
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    astre = 0.05;
    A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
    W0 = 4.0 * gamma0 / delta;           //ペナルティー項の係数[式(4.40)]
    M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]
    F0 = 50.0 / RR / temp;               //粒界移動の駆動力

    for (i = 1; i <= nm; i++)
    {
        for (j = 1; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            fij[i][j] = 0.0;
            anij[i][j] = 0;
            thij[i][j] = 0.0;
            if ((i == nm) || (j == nm))
            {
                fij[i][j] = F0;
                anij[i][j] = 1;
            }
            if (i > j)
            {
                fij[i][j] = -fij[i][j];
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                fij[i][j] = 0.0;
                anij[i][j] = 0;
            }
        }
    }
    // thij[1][3] = PI / 8.0;
    // thij[3][1] = PI / 8.0;

    double phi[N][NDX][NDY], phi2[N][NDX][NDY]; //フェーズフィールド、フェーズフィールド補助配列
    int phiIdx[N][NDX][NDY];                    //位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の番号
    int phiNum[NDX][NDY];

    for (k = 1; k <= nm; k++)
    {
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                phi[k][i][j] = 0.0;
                phi2[k][i][j] = 0.0; //フェーズフィールド、および補助配列の初期化
            }
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (ii = 1; ii <= nm - 1; ii++)
            {
                phi[ii][i][j] = 0.0;
            }
            phi[nm][i][j] = 1.0; // nm番目のフェーズフィールドを１に初期化
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            // if (i * i + j * j < 100)
            if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < 100)
            {
                phi[1][i][j] = 1.0;
                phi[2][i][j] = 0.0;
            }
            else
            {
                phi[1][i][j] = 0.0;
                phi[2][i][j] = 1.0;
            }
        }
    }
    // }

start:;

    if ((((int)(istep) % 100) == 0))
    {
        FILE *fp;

        fp = fopen("phi.dat", "w");

        for (int i = 0; i <= ndmx; i++)
        {
            for (int j = 0; j <= ndmy; j++)
            {
                fprintf(fp, "%e\n", phi[2][i][j]);
            }
        }
        fclose(fp);
    }

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

            //--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
            phinum = 0;
            for (ii = 1; ii <= nm; ii++)
            {
                if ((phi[ii][i][j] > 0.0) ||
                    ((phi[ii][i][j] == 0.0) && (phi[ii][ip][j] > 0.0) ||
                     (phi[ii][im][j] > 0.0) ||
                     (phi[ii][i][jp] > 0.0) ||
                     (phi[ii][i][jm] > 0.0)))
                {
                    phinum++;
                    phiIdx[phinum][i][j] = ii;
                    // printf("%d  ", n00);
                }
            }
            phiNum[i][j] = phinum;
        }
    }

    // Evolution Equations
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
                ip = 0;
            }
            if (i == 0)
            {
                im = ndmx;
            } //周期的境界条件
            if (j == ndmy)
            {
                jp = 0;
            }
            if (j == 0)
            {
                jm = ndmy;
            }

            for (n1 = 1; n1 <= phiNum[i][j]; n1++)
            {
                ii = phiIdx[n1][i][j];
                pddtt = 0.0;
                for (n2 = 1; n2 <= phiNum[i][j]; n2++)
                {
                    jj = phiIdx[n2][i][j];
                    sum1 = 0.0;
                    for (n3 = 1; n3 <= phiNum[i][j]; n3++)
                    {
                        kk = phiIdx[n3][i][j];

                        phidx = (phi[kk][ip][j] - phi[kk][im][j]) / 2.0; //フェーズフィールドの空間１階微分
                        phidy = (phi[kk][i][jp] - phi[kk][i][jm]) / 2.0;
                        phidxx = (phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j]); //フェーズフィールドの空間２階微分
                        phidyy = (phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j]);
                        phidxy = (phi[kk][ip][jp] + phi[kk][im][jm] - phi[kk][im][jp] - phi[kk][ip][jm]) / 4.0;

                        if (anij[ii][kk] == 1 && phidx * phidx + phidy * phidy != 0.0)
                        {
                            epsilon0 = sqrt(aij[ii][kk]);
                            theta0 = thij[ii][kk];

                            phidxp = phidx * cos(theta0) + phidy * sin(theta0);
                            phidyp = -phidx * sin(theta0) + phidy * cos(theta0);

                            phidxpx = phidxx * cos(theta0) + phidxy * sin(theta0);
                            phidypx = -phidxx * sin(theta0) + phidxy * cos(theta0);
                            phidxpy = phidxy * cos(theta0) + phidyy * sin(theta0);
                            phidypy = -phidxy * sin(theta0) + phidyy * cos(theta0);

                            phiabs = phidx * phidx + phidy * phidy;

                            ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 2.0));

                            epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy) * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0));
                            epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy) * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0));

                            term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep;

                            termx0 = (pow(phidxp, 3.0) * cos(theta0) - pow(phidyp, 3.0) * sin(theta0)) / phiabs;
                            termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0)) * phidx / pow(phiabs, 2.0);

                            termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * cos(theta0) - 3.0 * pow(phidyp, 2.0) * phidypx * sin(theta0)) / phiabs -
                                       (2.0 * phidx * phidxx + 2.0 * phidy * phidxy) * (pow(phidxp, 3.0) * cos(theta0) - pow(phidyp, 3.0) * sin(theta0)) / pow(phiabs, 2.0);

                            termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx))) / pow(phiabs, 2.0) -
                                       4.0 * (phidx * phidxx + phidy * phidxy) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0);

                            termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));

                            termy0 = (pow(phidxp, 3.0) * sin(theta0) + pow(phidyp, 3.0) * cos(theta0)) / phiabs;
                            termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0)) * phidy / pow(phiabs, 2.0);

                            termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * sin(theta0) + 3.0 * pow(phidyp, 2.0) * phidypy * cos(theta0)) / phiabs -
                                       (2.0 * phidx * phidxy + 2.0 * phidy * phidyy) * (pow(phidxp, 3.0) * sin(theta0) + pow(phidyp, 3.0) * cos(theta0)) / pow(phiabs, 2.0);

                            termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy))) / pow(phiabs, 2.0) -
                                       4.0 * (phidx * phidxy + phidy * phidyy) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0);

                            termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));

                            termiikk = term0 + termx + termy;
                        }
                        else
                        {
                            termiikk = aij[ii][kk] * (phidxx + phidyy);
                        }

                        if (anij[jj][kk] == 1 && phidx * phidx + phidy * phidy != 0.0)
                        {
                            epsilon0 = sqrt(aij[jj][kk]);
                            theta0 = thij[jj][kk];
                            phidxp = phidx * cos(theta0) + phidy * sin(theta0);
                            phidyp = -phidx * sin(theta0) + phidy * cos(theta0);

                            phidxpx = phidxx * cos(theta0) + phidxy * sin(theta0);
                            phidypx = -phidxx * sin(theta0) + phidxy * cos(theta0);
                            phidxpy = phidxy * cos(theta0) + phidyy * sin(theta0);
                            phidypy = -phidxy * sin(theta0) + phidyy * cos(theta0);

                            phiabs = phidx * phidx + phidy * phidy;

                            ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 2.0));

                            epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy) * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0));
                            epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy) * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0));

                            term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep;

                            termx0 = (pow(phidxp, 3.0) * cos(theta0) - pow(phidyp, 3.0) * sin(theta0)) / phiabs;
                            termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0)) * phidx / pow(phiabs, 2.0);

                            termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * cos(theta0) - 3.0 * pow(phidyp, 2.0) * phidypx * sin(theta0)) / phiabs -
                                       (2.0 * phidx * phidxx + 2.0 * phidy * phidxy) * (pow(phidxp, 3.0) * cos(theta0) - pow(phidyp, 3.0) * sin(theta0)) / pow(phiabs, 2.0);

                            termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx))) / pow(phiabs, 2.0) -
                                       4.0 * (phidx * phidxx + phidy * phidxy) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0);

                            termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));

                            termy0 = (pow(phidxp, 3.0) * sin(theta0) + pow(phidyp, 3.0) * cos(theta0)) / phiabs;
                            termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0)) * phidy / pow(phiabs, 2.0);

                            termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * sin(theta0) + 3.0 * pow(phidyp, 2.0) * phidypy * cos(theta0)) / phiabs -
                                       (2.0 * phidx * phidxy + 2.0 * phidy * phidyy) * (pow(phidxp, 3.0) * sin(theta0) + pow(phidyp, 3.0) * cos(theta0)) / pow(phiabs, 2.0);

                            termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy))) / pow(phiabs, 2.0) -
                                       4.0 * (phidx * phidxy + phidy * phidyy) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0)) / pow(phiabs, 3.0);

                            termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));

                            termjjkk = term0 + termx + termy;
                        }
                        else
                        {
                            termjjkk = aij[jj][kk] * (phidxx + phidyy);
                        }

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                    }
                    pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j] * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt(phi[ii][i][j] * phi[jj][i][j]));
                    //フェーズフィールドの発展方程式[式(4.31)]
                }
                phi2[ii][i][j] = phi[ii][i][j] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
                if (phi2[ii][i][j] >= 1.0)
                {
                    phi2[ii][i][j] = 1.0;
                } //フェーズフィールドの変域補正
                if (phi2[ii][i][j] <= 0.0)
                {
                    phi2[ii][i][j] = 0.0;
                }
            }
        } // j
    }     // i

    for (k = 1; k <= nm; k++)
    {
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                phi[k][i][j] = phi2[k][i][j];
            }
        }
    }

    //
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            sum1 = 0.0;
            for (k = 1; k <= nm; k++)
            {
                sum1 += phi[k][i][j];
            }
            for (k = 1; k <= nm; k++)
            {
                phi[k][i][j] = phi[k][i][j] / sum1;
            }
        }
    }

    istep = istep + 1.;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}
