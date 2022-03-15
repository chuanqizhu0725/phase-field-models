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

#define NDX 30 //差分計算における計算領域一辺の分割数
#define NDY 30 //差分計算における計算領域一辺の分割数
#define NDZ 30
#define N 3 //考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1; //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int ndmz = NDZ - 1;
int nm = N - 1, nmm = N - 2; //考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
double PI = 3.141592;        //π、計算カウント数
double RR = 8.3145;          //ガス定数

double aij[N][N]; //勾配エネルギー係数
double wij[N][N]; //ペナルティー項の係数
double mij[N][N]; //粒界の易動度
double fij[N][N]; //粒界移動の駆動力
int anij[N][N];
double thij[N][N];
double vphij[N][N];
int phinum;

int i, j, k, l, ii, jj, kk, ll, it; //整数
int ip, im, jp, jm, lp, lm;         //整数
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
double theta, thetadx, thetady, thetadz;
double vphi, vphidx, vphidy, vphidz;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidxx, phidxy, phidxz;
double phidy, phidyy, phidyz;
double phidz, phidzz;
double ep, ep1p, ep2p;
double epv, epvv, epvt, ept, eptt, eptv;

double term0, termx0, termx1, termy0, termy1, termz;

int x11, y11, x1h[10], y1h[10]; //初期核の座標
double t, r0, r;

double calcTheta(double dy, double dx);
double calcVphi(double dx, double dy, double dz);

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    nstep = 1001;
    dtime = 5.0;
    temp = 1000.0;
    L = 2000.0;
    vm0 = 7.0e-6;
    delta = 7.0;
    mobi = 1.0;

    dx = L / (double)NDX * 1.0e-9;       //差分プロック１辺の長さ(m)
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    astre = 0.1;
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
            vphij[i][j] = 0.0;
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

    double phi[N][NDX][NDY][NDZ], phi2[N][NDX][NDY][NDZ]; //フェーズフィールド、フェーズフィールド補助配列
    int phiIdx[N][NDX][NDY][NDZ];                         //位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の番号
    int phiNum[NDX][NDY][NDZ];

    for (k = 1; k <= nm; k++)
    {
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (l = 0; l <= ndmz; l++)
                {
                    phi[k][i][j][l] = 0.0;
                    phi2[k][i][j][l] = 0.0; //フェーズフィールド、および補助配列の初期化
                }
            }
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {
                for (ii = 1; ii <= nm - 1; ii++)
                {
                    phi[ii][i][j][l] = 0.0;
                }
                phi[nm][i][j][l] = 1.0; // nm番目のフェーズフィールドを１に初期化
            }
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {

                if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < 25)
                {

                    phi[1][i][j][l] = 1.0;
                    phi[2][i][j][l] = 0.0;
                }
                else
                {
                    phi[1][i][j][l] = 0.0;
                    phi[2][i][j][l] = 1.0;
                }
            }
        }
    }
    // }

start:;

    if ((((int)(istep) % 100) == 0))
    {
        FILE *stream;
        char buffer[30];
        sprintf(buffer, "3d%d.vtk", istep);
        stream = fopen(buffer, "a");

        fprintf(stream, "# vtk DataFile Version 1.0\n");
        fprintf(stream, "phi_%e.vtk\n", istep);
        fprintf(stream, "ASCII\n");
        fprintf(stream, "DATASET STRUCTURED_POINTS\n");
        fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
        fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
        fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
        fprintf(stream, "\n");
        fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
        fprintf(stream, "SCALARS scalars float\n");
        fprintf(stream, "LOOKUP_TABLE default\n");

        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (l = 0; l <= ndmz; l++)
                {
                    sum2 = 0.0;
                    for (k = 0; k <= nm; k++)
                    {
                        sum2 += phi[k][i][j][l] * phi[k][i][j][l];
                    }
                    fprintf(stream, "%e\n", sum2);
                }
            }
        }
        fclose(stream);
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                lp = l + 1;
                lm = l - 1;
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
                if (l == ndmz)
                {
                    lp = 0;
                }
                if (l == 0)
                {
                    lm = ndmz;
                }

                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if ((phi[ii][i][j][l] > 0.0) ||
                        ((phi[ii][i][j][l] == 0.0) && (phi[ii][ip][j][l] > 0.0) ||
                         (phi[ii][im][j][l] > 0.0) ||
                         (phi[ii][i][jp][l] > 0.0) ||
                         (phi[ii][i][jm][l] > 0.0) ||
                         (phi[ii][i][j][lp] > 0.0) ||
                         (phi[ii][i][j][lm] > 0.0)))
                    {
                        phinum++;
                        phiIdx[phinum][i][j][l] = ii;
                    }
                }
                phiNum[i][j][l] = phinum;
            }
        }
    }

    // Evolution Equations
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                lp = l + 1;
                lm = l - 1;
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
                if (l == ndmz)
                {
                    lp = 0;
                }
                if (l == 0)
                {
                    lm = ndmz;
                }

                for (n1 = 1; n1 <= phiNum[i][j][l]; n1++)
                {
                    ii = phiIdx[n1][i][j][l];
                    pddtt = 0.0;
                    for (n2 = 1; n2 <= phiNum[i][j][l]; n2++)
                    {
                        jj = phiIdx[n2][i][j][l];
                        sum1 = 0.0;
                        for (n3 = 1; n3 <= phiNum[i][j][l]; n3++)
                        {
                            kk = phiIdx[n3][i][j][l];

                            phidx = (phi[kk][ip][j][l] - phi[kk][im][j][l]) / 2.0; //フェーズフィールドの空間１階微分
                            phidy = (phi[kk][i][jp][l] - phi[kk][i][jm][l]) / 2.0;
                            phidz = (phi[kk][i][j][lp] - phi[kk][i][j][lm]) / 2.0;

                            phidxx = (phi[kk][ip][j][l] + phi[kk][im][j][l] - 2.0 * phi[kk][i][j][l]); //フェーズフィールドの空間２階微分
                            phidyy = (phi[kk][i][jp][l] + phi[kk][i][jm][l] - 2.0 * phi[kk][i][j][l]);
                            phidzz = (phi[kk][i][j][lp] + phi[kk][i][j][lm] - 2.0 * phi[kk][i][j][l]);

                            phidxy = (phi[kk][ip][jp][l] + phi[kk][im][jm][l] - phi[kk][im][jp][l] - phi[kk][ip][jm][l]) / 4.0;
                            phidxz = (phi[kk][ip][j][lp] + phi[kk][im][j][lm] - phi[kk][im][j][lp] - phi[kk][ip][j][lm]) / 4.0;
                            phidyz = (phi[kk][i][jp][lp] + phi[kk][i][jm][lm] - phi[kk][i][jm][lp] - phi[kk][i][jp][lm]) / 4.0;

                            if (anij[ii][kk] == 1)
                            {
                                if (phidx * phidx + phidy * phidy != 0.0)
                                {
                                    theta = calcTheta(phidy, phidx) - thij[ii][kk];
                                    thetadx = (phidxy * phidx - phidxx * phidy) / (phidx * phidx + phidy * phidy);
                                    thetady = (phidyy * phidx - phidxy * phidy) / (phidx * phidx + phidy * phidy);
                                    thetadz = (phidyz * phidx - phidxz * phidy) / (phidx * phidx + phidy * phidy);

                                    // vphi = calcVphi(phidx, phidy, phidz) - vphij[ii][kk];
                                    vphi = PI / 2.0;
                                    vphidx = ((phidxx * cos(theta) + phidxy * sin(theta)) * phidz + sqrt(phidx * phidx + phidy * phidy) * phidxz) / (phidx * phidx + phidy * phidy + phidz * phidz);
                                    vphidy = ((phidxy * cos(theta) + phidyy * sin(theta)) * phidz + sqrt(phidx * phidx + phidy * phidy) * phidyz) / (phidx * phidx + phidy * phidy + phidz * phidz);
                                    vphidz = ((phidxz * cos(theta) + phidyz * sin(theta)) * phidz + sqrt(phidx * phidx + phidy * phidy) * phidzz) / (phidx * phidx + phidy * phidy + phidz * phidz);

                                    epsilon0 = sqrt(aij[ii][kk]);
                                    ep = epsilon0 * (1 - 3.0 * astre + 4.0 * astre * (pow(sin(vphi), 4.0) * (pow(cos(theta), 4.0) + pow(sin(theta), 4.0)) + pow(cos(vphi), 4.0)));
                                    ept = -4.0 * epsilon0 * astre * pow(sin(vphi), 4.0) * sin(4.0 * theta);
                                    eptt = -16.0 * epsilon0 * astre * pow(sin(vphi), 4.0) * cos(4.0 * theta);
                                    eptv = -16.0 * epsilon0 * astre * pow(sin(vphi), 3.0) * sin(4.0 * theta) * cos(vphi);
                                    epv = 16.0 * epsilon0 * astre * (pow(sin(vphi), 3.0) * cos(vphi) * (pow(cos(theta), 4.0) + pow(sin(theta), 4.0)) - pow(cos(vphi), 3.0) * sin(vphi));
                                    epvv = 16.0 * epsilon0 * astre * (3.0 * pow(sin(vphi), 2.0) * pow(cos(vphi), 2.0) * (pow(cos(theta), 4.0) + pow(sin(theta), 4.0) + 1) - pow(sin(vphi), 4.0) * (pow(sin(theta), 4.0) + pow(cos(theta), 4.0)) - pow(cos(vphi), 4.0));
                                    epvt = 16.0 * epsilon0 * astre * pow(sin(vphi), 3.0) * cos(vphi) * sin(4.0 * theta);

                                    term0 = ep * ep * (phidxx + phidyy + phidzz) + 2.0 * ep * ept * (phidx * thetadx + phidy * thetady + phidz * thetadz) + 2.0 * ep * epv * (phidx * vphidx + phidy * vphidy + phidz * vphidz);

                                    termx0 = (((epv * vphidx + ept * thetadx) * epv + (epvv * vphidx + epvt * thetadx) * ep) * cos(theta) - sin(theta) * thetadx * ep * epv) * phidz + phidxz * ep * epv * cos(theta);
                                    termy0 = (((epv * vphidy + ept * thetady) * epv + (epvv * vphidy + epvt * thetady) * ep) * sin(theta) + cos(theta) * thetady * ep * epv) * phidz + phidyz * ep * epv * sin(theta);

                                    termx1 = -(((epv * vphidx + ept * thetadx) * ept + (eptv * vphidx + eptt * thetadx) * ep) * phidy + phidxy * ep * ept) / pow(sin(vphi), 2.0) + 2.0 * cos(vphi) / pow(sin(vphi), 3.0) * vphidx * ep * ept * phidy;
                                    termy1 = (((epv * vphidy + ept * thetady) * ept + (eptv * vphidy + eptt * thetady) * ep) * phidx + phidxy * ep * ept) / pow(sin(vphi), 2.0) - 2.0 * cos(vphi) / pow(sin(vphi), 3.0) * vphidy * ep * ept * phidx;

                                    termz = -((epv * vphidz + ept * thetadz) * epv + ep * (epvv * vphidz + epvt * thetadz)) * sqrt(phidx * phidx + phidy * phidy) - 1.0 / sqrt(phidx * phidx + phidy * phidy) * (phidx * phidxz + phidy * phidyz) * ep * epv;

                                    termiikk = term0 + termx0 + termx1 + termy0 + termy1 + termz;
                                }
                                else
                                {
                                    termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                                }
                            }
                            else
                            {
                                termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                            }

                            if (anij[jj][kk] == 1)
                            {
                                if (phidx * phidx + phidy * phidy != 0.0)
                                {
                                    theta = calcTheta(phidy, phidx) - thij[ii][kk];
                                    thetadx = (phidxy * phidx - phidxx * phidy) / (phidx * phidx + phidy * phidy);
                                    thetady = (phidyy * phidx - phidxy * phidy) / (phidx * phidx + phidy * phidy);
                                    thetadz = (phidyz * phidx - phidxz * phidy) / (phidx * phidx + phidy * phidy);

                                    // vphi = calcVphi(phidx, phidy, phidz) - vphij[ii][kk];
                                    vphi = PI / 2.0;
                                    vphidx = ((phidxx * cos(theta) + phidxy * sin(theta)) * phidz + sqrt(phidx * phidx + phidy * phidy) * phidxz) / (phidx * phidx + phidy * phidy + phidz * phidz);
                                    vphidy = ((phidxy * cos(theta) + phidyy * sin(theta)) * phidz + sqrt(phidx * phidx + phidy * phidy) * phidyz) / (phidx * phidx + phidy * phidy + phidz * phidz);
                                    vphidz = ((phidxz * cos(theta) + phidyz * sin(theta)) * phidz + sqrt(phidx * phidx + phidy * phidy) * phidzz) / (phidx * phidx + phidy * phidy + phidz * phidz);

                                    epsilon0 = sqrt(aij[jj][kk]);
                                    ep = epsilon0 * (1 - 3.0 * astre + 4.0 * astre * (pow(sin(vphi), 4.0) * (pow(cos(theta), 4.0) + pow(sin(theta), 4.0)) + pow(cos(vphi), 4.0)));
                                    ept = -4.0 * epsilon0 * astre * pow(sin(vphi), 4.0) * sin(4.0 * theta);
                                    eptt = -16.0 * epsilon0 * astre * pow(sin(vphi), 4.0) * cos(4.0 * theta);
                                    eptv = -16.0 * epsilon0 * astre * pow(sin(vphi), 3.0) * sin(4.0 * theta);
                                    epv = 16.0 * epsilon0 * astre * (pow(sin(vphi), 3.0) * cos(vphi) * (pow(cos(theta), 4.0) + pow(sin(theta), 4.0)) - pow(cos(vphi), 3.0) * sin(vphi));
                                    epvv = 16.0 * epsilon0 * astre * (3.0 * pow(sin(vphi), 2.0) * pow(cos(vphi), 2.0) * (pow(cos(theta), 4.0) + pow(sin(theta), 4.0) + 1) - pow(sin(vphi), 4.0) * (pow(sin(theta), 4.0) + pow(cos(theta), 4.0)) - pow(cos(vphi), 4.0));
                                    epvt = 16.0 * epsilon0 * astre * pow(sin(vphi), 3.0) * cos(vphi) * sin(4.0 * theta);

                                    term0 = ep * ep * (phidxx + phidyy + phidzz) + 2.0 * ep * ept * (phidx * thetadx + phidy * thetady + phidz * thetadz) + 2.0 * ep * epv * (phidx * vphidx + phidy * vphidy + phidz * vphidz);
                                    termx0 = (((epv * vphidx + ept * thetadx) * epv + (epvv * vphidx + epvt * thetadx) * ep) * cos(theta) - sin(theta) * thetadx * ep * epv) * phidz + phidxz * ep * epv * cos(theta);
                                    termx1 = -(((epv * vphidx + ept * thetadx) * ept + (eptv * vphidx + eptt * thetadx) * ep) * phidy + phidxy * ep * ept) / pow(sin(vphi), 2.0) + 2.0 * cos(vphi) / pow(sin(vphi), 3.0) * vphidx * ep * ept * phidy;
                                    termy0 = (((epv * vphidy + ept * thetady) * epv + (epvv * vphidy + epvt * thetady) * ep) * sin(theta) + cos(theta) * thetady * ep * epv) * phidz + phidyz * ep * epv * sin(theta);
                                    termy1 = (((epv * vphidy + ept * thetady) * ept + (eptv * vphidy + eptt * thetady) * ep) * phidx + phidxy * ep * ept) / pow(sin(vphi), 2.0) - 2.0 * cos(vphi) / pow(sin(vphi), 3.0) * vphidy * ep * ept * phidx;
                                    termz = -((epv * vphidz + ept * thetadz) * epv + ep * (epvv * vphidz + epvt * thetadz)) * sqrt(phidx * phidx + phidy * phidy) - 1.0 / sqrt(phidx * phidx + phidy * phidy) * (phidx * phidxz + phidy * phidyz) * ep * epv;

                                    termjjkk = term0 + termx0 + termx1 + termy0 + termy1 + termz;
                                }
                                else
                                {
                                    termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                                }
                            }
                            else
                            {
                                termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                            }

                            sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][l];
                        }
                        pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j][l] * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt(phi[ii][i][j][l] * phi[jj][i][j][l]));
                        //フェーズフィールドの発展方程式[式(4.31)]
                    }
                    phi2[ii][i][j][l] = phi[ii][i][j][l] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
                    if (phi2[ii][i][j][l] >= 1.0)
                    {
                        phi2[ii][i][j][l] = 1.0;
                    } //フェーズフィールドの変域補正
                    if (phi2[ii][i][j][l] <= 0.0)
                    {
                        phi2[ii][i][j][l] = 0.0;
                    }
                }
            } // j
        }     // i
    }

    for (k = 1; k <= nm; k++)
    {
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (l = 0; l <= ndmz; l++)
                {
                    phi[k][i][j][l] = phi2[k][i][j][l];
                }
            }
        }
    }

    //
    // for (i = 0; i <= ndmx; i++)
    // {
    //     for (j = 0; j <= ndmy; j++)
    //     {
    //         for (l = 0; l <= ndmz; l++)
    //         {
    //             sum1 = 0.0;
    //             for (k = 1; k <= nm; k++)
    //             {
    //                 sum1 += phi[k][i][j][l];
    //             }
    //             for (k = 1; k <= nm; k++)
    //             {
    //                 phi[k][i][j][l] = phi[k][i][j][l] / sum1;
    //             }
    //         }
    //     }
    // }

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

double calcVphi(double dx, double dy, double dz)
{
    if (dz != 0.0)
    {
        return atan(sqrt(dx * dx + dy * dy) / dz);
    }
    else if (dz == 0.0)
    {
        return PI / 2.0;
    }
    else
    {
        return 0;
    }
}
