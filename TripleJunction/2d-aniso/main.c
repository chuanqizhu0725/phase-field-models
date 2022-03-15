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

#define NDX 100 //差分計算における計算領域一辺の分割数
#define NDY 100 //差分計算における計算領域一辺の分割数
#define N 4     //考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)

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
double ep, ep1p, ep2p;

int x11, y11, x1h[10], y1h[10]; //初期核の座標
double t, r0, r;

double calcTheta(double dy, double dx);

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    nstep = 10001;
    dtime = 5.0;
    temp = 1000.0;
    L = 2000.0;
    vm0 = 7.0e-6;
    delta = 7.0;
    mobi = 1.0;

    dx = L / (double)NDX * 1.0e-9;       //差分プロック１辺の長さ(m)
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    astre = 0.03;
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
    thij[1][3] = PI / 8.0;
    thij[3][1] = PI / 8.0;

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
            if (i > NDX * 4 / 5 && j < NDY / 2)
            {
                phi[1][i][j] = 1.0;
                phi[2][i][j] = 0.0;
                phi[3][i][j] = 0.0;
            }
            else if (i > NDX * 4 / 5 && j >= NDY / 2)
            {
                phi[1][i][j] = 0.0;
                phi[2][i][j] = 1.0;
                phi[3][i][j] = 0.0;
            }
            else
            {
                phi[1][i][j] = 0.0;
                phi[2][i][j] = 0.0;
                phi[3][i][j] = 1.0;
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
                sum2 = 0.0;
                for (int k = 0; k <= nm; k++)
                {
                    sum2 += phi[k][i][j] * phi[k][i][j];
                }
                fprintf(fp, "%e\n", sum2);
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
                ip = ndmx;
            }
            if (i == 0)
            {
                im = 0;
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
                ip = ndmx;
            }
            if (i == 0)
            {
                im = 0;
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

                        if (anij[ii][kk] == 1)
                        {
                            epsilon0 = sqrt(aij[ii][kk]);
                            theta0 = thij[ii][kk];
                            theta = calcTheta(phidy, phidx);                                    //界面の法線方向の角度[式(4.24)]
                            ep = epsilon0 * (1.0 + astre * cos(4.0 * (theta - theta0)));        //勾配エネルギー係数の平方根[式(4.23)]
                            ep1p = -epsilon0 * astre * 4.0 * sin(4.0 * (theta - theta0));       // epの角度による１階微分
                            ep2p = -epsilon0 * astre * 4.0 * 4.0 * cos(4.0 * (theta - theta0)); // epの角度による２階微分

                            thetax = (phidxy * phidx - phidxx * phidy) / (phidx * phidx + phidy * phidy);
                            thetay = (phidyy * phidx - phidxy * phidy) / (phidx * phidx + phidy * phidy);

                            termiikk = ep * ep * (phidxx + phidyy) +
                                       2.0 * ep * ep1p * (phidx * thetax + phidy * thetay) -
                                       (ep1p * ep1p + ep * ep2p) * (phidxy * phidx * phidy - phidxx * phidy * phidy) / (phidx * phidx + phidy * phidy) +
                                       (ep1p * ep1p + ep * ep2p) * (phidyy * phidx * phidx - phidxy * phidx * phidy) / (phidx * phidx + phidy * phidy);

                            // termiikk = ep * ep * (phidxx + phidyy) + ep * ep1p * ((phidyy - phidxx) * sin(2.0 * theta) + 2.0 * phidxy * cos(2.0 * theta)) - 0.5 * (ep1p * ep1p + ep * ep2p) * (2.0 * phidxy * sin(2.0 * theta) - phidxx - phidyy - (phidyy - phidxx) * cos(2.0 * theta));
                        }
                        else
                        {
                            termiikk = aij[ii][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]);
                        }

                        if (anij[jj][kk] == 1)
                        {
                            epsilon0 = sqrt(aij[jj][kk]);
                            theta0 = thij[jj][kk];
                            theta = calcTheta(phidy, phidx);                                    //界面の法線方向の角度[式(4.24)]
                            ep = epsilon0 * (1.0 + astre * cos(4.0 * (theta - theta0)));        //勾配エネルギー係数の平方根[式(4.23)]
                            ep1p = -epsilon0 * astre * 4.0 * sin(4.0 * (theta - theta0));       // epの角度による１階微分
                            ep2p = -epsilon0 * astre * 4.0 * 4.0 * cos(4.0 * (theta - theta0)); // epの角度による２階微分

                            thetax = (phidxy * phidx - phidxx * phidy) / (phidx * phidx + phidy * phidy);
                            thetay = (phidyy * phidx - phidxy * phidy) / (phidx * phidx + phidy * phidy);

                            termjjkk = ep * ep * (phidxx + phidyy) +
                                       2.0 * ep * ep1p * (phidx * thetax + phidy * thetay) -
                                       (ep1p * ep1p + ep * ep2p) * (phidxy * phidx * phidy - phidxx * phidy * phidy) / (phidx * phidx + phidy * phidy) +
                                       (ep1p * ep1p + ep * ep2p) * (phidyy * phidx * phidx - phidxy * phidx * phidy) / (phidx * phidx + phidy * phidy);
                        }
                        else
                        {
                            termjjkk = aij[jj][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]);
                        }

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                    }
                    pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j]) * (sum1 - 8.0 / PI * fij[ii][jj] * (0.5 * i * i / NDX / NDX) * sqrt(phi[ii][i][j] * phi[jj][i][j]));
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

double calcTheta(double dy, double dx)
{
    if (dx != 0.0)
    {
        return atan(dy / dx);
    }
    else if (dx == 0.0 && dy > 0)
    {
        return PI / 2.0;
    }
    else if (dx == 0.0 && dy < 0)
    {
        return PI / (-2.0);
    }
    else
    {
        return 0;
    }
}
