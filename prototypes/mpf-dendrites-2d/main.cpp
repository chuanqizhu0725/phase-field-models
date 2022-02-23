// Adapted from polycrystal-binary
// Author: Chuanqi Zhu
// Created on: 2022/2/22

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

    A_alph = 1.0e+02 / RR / temp;
    c_alph0 = 0.1; //化学的自由エネルギー内のパラメータ
    A_beta = 1.0e+02 / RR / temp;
    c_beta0 = 0.4;
    Da = 0.001; //各相の拡散係数は全て等しいと仮定
    Db = 0.005;

    initialize();

    log();

start:;

    if ((((int)(istep) % pstep) == 0))
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

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            c = ch[i][j]; //濃度場
            sum1 = 0.0;
            for (ii = 1; ii <= nm - 1; ii++)
            {
                sum1 += phi[ii][i][j];
            }
            sa = sah[i][j] = sum1; // Sαの計算[式(4.46)]
            sum1 = 0.0;
            for (ii = nm; ii <= nm; ii++)
            {
                sum1 += phi[ii][i][j];
            }
            sb = sbh[i][j] = sum1; // Sβの計算[式(4.46)]

            //局所平衡組成の計算[式(4.47)]
            cah[i][j] = (A_beta * c + (A_alph * c_alph0 - A_beta * c_beta0) * sb) / (A_beta * sa + A_alph * sb);
            cbh[i][j] = (A_alph * c + (A_beta * c_beta0 - A_alph * c_alph0) * sa) / (A_beta * sa + A_alph * sb);
            if (cah[i][j] >= 1.0)
            {
                cah[i][j] = 1.0;
            }
            if (cah[i][j] <= 0.0)
            {
                cah[i][j] = 0.0;
            } //濃度場の変域補正
            if (cbh[i][j] >= 1.0)
            {
                cbh[i][j] = 1.0;
            }
            if (cbh[i][j] <= 0.0)
            {
                cbh[i][j] = 0.0;
            }
        }
    }

    // Evolution Equation of Phase Fields
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
                    if (ii != nm && jj == nm)
                    {
                        // cii = cah[i][j];
                        // cjj = cbh[i][j];
                        // gii = A_alph * (cii - c_alph0) * (cii - c_alph0);
                        // gjj = A_beta * (cjj - c_beta0) * (cjj - c_beta0) + 40.0 / RR / temp;
                        // ptl = 2.0 * A_alph * (cii - c_alph0);
                        // dF = gii - gjj - ptl * (cii - cjj);
                        dF = (450.0 - 1200.0 * cbh[i][j]) / RR / temp;
                    }
                    else if (ii == nm && jj != nm)
                    {
                        // cii = cbh[i][j];
                        // cjj = cah[i][j];
                        // gii = A_beta * (cii - c_beta0) * (cii - c_beta0) + 40.0 / RR / temp;
                        // gjj = A_alph * (cjj - c_alph0) * (cjj - c_alph0);
                        // ptl = 2.0 * A_alph * (cjj - c_alph0);
                        // dF = gii - gjj - ptl * (cii - cjj);
                        dF = -(450.0 - 1200.0 * cbh[i][j]) / RR / temp;
                    }
                    else
                    {
                        dF = 0.0;
                    }
                    // double dF = gii - gjj;
                    pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j] * phi[jj][i][j]));
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

    // Evolution Equation of Concentration field
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

            //拡散方程式内における微分計算
            dev1_a = 0.25 * ((sah[ip][j] - sah[im][j]) * (cah[ip][j] - cah[im][j]) + (sah[i][jp] - sah[i][jm]) * (cah[i][jp] - cah[i][jm]));
            dev1_b = 0.25 * ((sbh[ip][j] - sbh[im][j]) * (cbh[ip][j] - cbh[im][j]) + (sbh[i][jp] - sbh[i][jm]) * (cbh[i][jp] - cbh[i][jm]));
            dev2_a = sah[i][j] * (cah[ip][j] + cah[im][j] + cah[i][jp] + cah[i][jm] - 4.0 * cah[i][j]);
            dev2_b = sbh[i][j] * (cbh[ip][j] + cbh[im][j] + cbh[i][jp] + cbh[i][jm] - 4.0 * cbh[i][j]);

            cddtt = Da * (dev1_a + dev2_a) + Db * (dev1_b + dev2_b); //拡散方程式[式(4.42)]
            // ch2[i][j]=ch[i][j]+cddtt*delt;	//濃度場の時間発展(陽解法)
            ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 1; k <= nm; k++)
            {
                phi[k][i][j] = phi2[k][i][j]; //補助配列を主配列に移動（フェーズフィールド）
            }
            ch[i][j] = ch2[i][j]; //補助配列を主配列に移動（濃度場）
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

    //*** 濃度場の収支補正 *************************************************************
    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            sum1 += ch[i][j];
        }
    }
    dc0 = sum1 / nd / nd - c0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            ch[i][j] = ch[i][j] - dc0;
            if (ch[i][j] > 1.0)
            {
                ch[i][j] = 1.0;
            }
            if (ch[i][j] < 0.0)
            {
                ch[i][j] = 0.0;
            }
        }
    }

    istep = istep + 1;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}
