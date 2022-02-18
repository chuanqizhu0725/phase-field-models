// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/2/16

#include "header.h"

//******* メインプログラム ******************************************
int main(void)
{
    datain();

    dx = L / (double)ND * 1.0e-9;        //差分プロック１辺の長さ(m)
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
    W0 = 4.0 * gamma0 / delta;           //ペナルティー項の係数[式(4.40)]
    M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]
    F0 = 80.0 / RR / temp;               //粒界移動の駆動力

    initialize();

    log();

start:;

    if ((((int)(istep) % 1000) == 0))
    {
        datasave(istep);
    }

    //**** 各差分プロックにおけるphiNum[i][j]とphiIdx[n00][i][j]を調査 *********************
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            ip = i + 1;
            im = i - 1;
            jp = j + 1;
            jm = j - 1;
            if (i == ndm)
            {
                ip = 0;
            }
            if (i == 0)
            {
                im = ndm;
            } //周期的境界条件
            if (j == ndm)
            {
                jp = 0;
            }
            if (j == 0)
            {
                jm = ndm;
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
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            ip = i + 1;
            ipp = i + 2;
            im = i - 1;
            imm = i - 2;
            jp = j + 1;
            jpp = j + 2;
            jm = j - 1;
            jmm = j - 2;
            if (i == ndm)
            {
                ip = 0;
                ipp = 1;
            }
            if (i == ndm - 1)
            {
                ipp = 0;
            }
            if (i == 0)
            {
                im = ndm;
                imm = ndm - 1;
            } //周期的境界条件
            if (i == 1)
            {
                imm = ndm;
            }
            if (j == ndm)
            {
                jp = 0;
                jpp = 1;
            }
            if (j == ndm - 1)
            {
                jpp = 0;
            }
            if (j == 0)
            {
                jm = ndm;
                jmm = ndm - 1;
            }
            if (j == 1)
            {
                jmm = ndm;
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
                        phidxx = phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j]; //フェーズフィールドの空間２階微分
                        phidyy = phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j];
                        phidxy = (phi[kk][ip][jp] + phi[kk][im][jm] - phi[kk][im][jp] - phi[kk][ip][jm]) / 4.0;

                        if (anij[ii][kk])
                        {
                            epsilon0 = sqrt(aij[ii][kk]);
                            theta0 = thij[ii][kk];
                            theta = calcTheta(phidy, phidx);                                    //界面の法線方向の角度[式(4.24)]
                            ep = epsilon0 * (1.0 + astre * cos(4.0 * (theta - theta0)));        //勾配エネルギー係数の平方根[式(4.23)]
                            ep1p = -epsilon0 * astre * 4.0 * sin(4.0 * (theta - theta0));       // epの角度による１階微分
                            ep2p = -epsilon0 * astre * 4.0 * 4.0 * cos(4.0 * (theta - theta0)); // epの角度による２階微分

                            termiikk = ep * ep * (phidxx + phidyy) + ep * ep1p * ((phidyy - phidxx) * sin(2.0 * theta) + 2.0 * phidxy * cos(2.0 * theta)) - 0.5 * (ep1p * ep1p + ep * ep2p) * (2.0 * phidxy * sin(2.0 * theta) - phidxx - phidyy - (phidyy - phidxx) * cos(2.0 * theta));
                        }
                        else
                        {
                            termiikk = aij[ii][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]);
                        }

                        if (anij[jj][kk])
                        {
                            epsilon0 = sqrt(aij[jj][kk]);
                            theta0 = thij[jj][kk];
                            theta = calcTheta(phidy, phidx);                                    //界面の法線方向の角度[式(4.24)]
                            ep = epsilon0 * (1.0 + astre * cos(4.0 * (theta - theta0)));        //勾配エネルギー係数の平方根[式(4.23)]
                            ep1p = -epsilon0 * astre * 4.0 * sin(4.0 * (theta - theta0));       // epの角度による１階微分
                            ep2p = -epsilon0 * astre * 4.0 * 4.0 * cos(4.0 * (theta - theta0)); // epの角度による２階微分

                            termjjkk = ep * ep * (phidxx + phidyy) + ep * ep1p * ((phidyy - phidxx) * sin(2.0 * theta) + 2.0 * phidxy * cos(2.0 * theta)) - 0.5 * (ep1p * ep1p + ep * ep2p) * (2.0 * phidxy * sin(2.0 * theta) - phidxx - phidyy - (phidyy - phidxx) * cos(2.0 * theta));
                        }
                        else
                        {
                            termjjkk = aij[jj][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]);
                        }

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j]; //[式(4.31)の一部]
                    }
                    pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j]) * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt(phi[ii][i][j] * phi[jj][i][j]));
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
        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                phi[k][i][j] = phi2[k][i][j];
            }
        }
    }

    //
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
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
