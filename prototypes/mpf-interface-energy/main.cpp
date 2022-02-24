#include "header.h"

int main(void)
{
    // datain();

    nstep = 201;
    pstep = 10;

    dx = 30.0e-9;   // m
    dtime = 4.0e-7; // s

    gamma0 = 0.2; // J/m2
    delta = 7.0 * dx;

    temp = 1000;                  // K
    vm0 = 7.0e-6;                 // mole/m3
    mobi = 5.0 * vm0 / RR / temp; // 4.20951e-09

    A0 = 8.0 * delta * gamma0 / PI / PI;
    W0 = 4.0 * gamma0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);

    cout << "The phase field mobility is: " << mobi << endl;
    cout << "The Gradient coefficient is: " << A0 << endl;
    cout << "The penalty coefficient is: " << W0 << endl;
    cout << "The stablity number is: " << dtime * M0 * A0 / (dx * dx) << endl;

    initialize();

start:;

    // calculate the excessive enenrgy at the interface region
    Eexc = 0.0;
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
            phidx = (phi[1][ip][j] - phi[1][im][j]) / (2.0 * dx);
            phidy = (phi[1][i][jp] - phi[1][i][jm]) / (2.0 * dx);
            Eexc += 0.5 * A0 * (phidx * phidx + phidy * phidy) + W0 * phi[1][i][j] * (1.0 - phi[1][i][j]);
        }
    }

    if ((((int)(istep) % pstep) == 0))
    {
        // datasave(istep);
        cout << "the excessive energy density (J/m2) is: " << 2.0 * Eexc * 7.0 * dx << endl;
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

                        termiikk = aij[ii][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]) / (dx * dx);

                        termjjkk = aij[jj][kk] * (phi[kk][ip][j] + phi[kk][im][j] + phi[kk][i][jp] + phi[kk][i][jm] - 4.0 * phi[kk][i][j]) / (dx * dx);

                        sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j];
                    }
                    pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j]) * sum1;
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

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 1; k <= nm; k++)
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

    istep = istep + 1;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}