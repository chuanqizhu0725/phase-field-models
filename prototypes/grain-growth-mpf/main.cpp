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
    F0 = 50.0 / RR / temp;               //粒界移動の駆動力

    initialize();

    log();
    clock_t start_t, end_t, total_t;
    start_t = clock();

start:;

    if ((((int)(istep) % 10000) == 0))
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
                        sum1 += 0.5 * (aij[ii][kk] - aij[jj][kk]) * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j]; //[式(4.31)の一部]
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
    end_t = clock();
    total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    printf("Total time taken: %lu secs\n", total_t);
    printf("Exiting of the program...\n");
    return 0;
}
