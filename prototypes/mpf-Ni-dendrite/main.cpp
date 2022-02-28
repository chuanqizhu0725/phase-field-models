#include "header.h"

int main(void)
{
    // datain();

    nstep = 100001;
    pstep = 1000;

    // dx = 30.0e-9;    // m
    // dtime = 1.0e-11; // s

    dx = 5.0e-7;    // m
    dtime = 3.0e-9; // s

    delta = 4.0 * dx;
    gamma0 = 0.37; // J/m2
    astre = 0.05;

    Tm = 1728.0; // K
    Tini = 1511.2;

    cndct = 84.01;    // thermal conductivity
    speht = 5.42e+06; // J/(K*m3)
    laht = 2.350e+09; // J/m3

    mobi = 4.20951e-05;

    A0 = 8.0 * delta * gamma0 / PI / PI;
    W0 = 4.0 * gamma0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);

    cout << "The stablity number of diffusion equaiton is: " << dtime / dx / dx * cndct / speht << endl;
    cout << "The stablity number of pf equaiton is: " << dtime / dx / dx * mobi * gamma0 << endl;
    cout << "The upper limit of driving force is: " << 4.0 * dx / dtime / mobi / PI << endl;

    initialize();

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
                ip = ndm - 1;
            }
            if (i == 0)
            {
                im = 1;
            }
            if (j == ndm)
            {
                jp = ndm - 1;
            }
            if (j == 0)
            {
                jm = 1;
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
                }
            }
            phiNum[i][j] = phinum;
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
                ip = ndm - 1;
            }
            if (i == 0)
            {
                im = 1;
            }
            if (j == ndm)
            {
                jp = ndm - 1;
            }
            if (j == 0)
            {
                jm = 1;
            }
            tddtt = 0.0;
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

                        phidx = (phi[kk][ip][j] - phi[kk][im][j]) / 2.0 / dx; //フェーズフィールドの空間１階微分
                        phidy = (phi[kk][i][jp] - phi[kk][i][jm]) / 2.0 / dx;
                        phidxx = (phi[kk][ip][j] + phi[kk][im][j] - 2.0 * phi[kk][i][j]) / (dx * dx); //フェーズフィールドの空間２階微分
                        phidyy = (phi[kk][i][jp] + phi[kk][i][jm] - 2.0 * phi[kk][i][j]) / (dx * dx);
                        phidxy = (phi[kk][ip][jp] + phi[kk][im][jm] - phi[kk][im][jp] - phi[kk][ip][jm]) / 4.0 / (dx * dx);

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
                            termiikk = aij[ii][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]) / (dx * dx);
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
                            termjjkk = aij[jj][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]) / (dx * dx);
                        }

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                    }
                    if (ii != nm && jj == nm)
                    {
                        // dF = 3.0e6;

                        dF = (Tm - temp[i][j]) / Tm * 2.0e6;
                    }
                    else if (ii == nm && jj != nm)
                    {
                        // dF = -3.0e6;
                        dF = -(Tm - temp[i][j]) / Tm * 2.0e6;
                    }
                    else
                    {
                        dF = 0.0;
                    }
                    double pdtiijj = -2.0 * mij[ii][jj] / double(phiNum[i][j]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j] * phi[jj][i][j]));
                    // release heat when solid phase grows into liquid
                    if (ii != nm && jj == nm)
                    {
                        tddtt += 30.0 * phi[ii][i][j] * (1.0 - phi[ii][i][j]) * phi[ii][i][j] * (1.0 - phi[ii][i][j]) * laht * pdtiijj / speht;
                    }
                    // else if (ii == nm && jj != nm) // release heat when liquid phase shrinks
                    // {
                    //     tddtt -= 30.0 * phi[ii][i][j] * (1.0 - phi[ii][i][j]) * phi[ii][i][j] * (1.0 - phi[ii][i][j]) * laht * pdtiijj / speht;
                    // }
                    pddtt += pdtiijj;
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

            tddtt += cndct / speht * (temp[ip][j] + temp[im][j] + temp[i][jp] + temp[i][jm] - 4.0 * temp[i][j]) / dx / dx;
            temp2[i][j] = temp[i][j] + tddtt * dtime;
        } // j
    }     // i

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 1; k <= nm; k++)
            {
                phi[k][i][j] = phi2[k][i][j];
            }
            temp[i][j] = temp2[i][j];
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

    istep = istep + 1;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}