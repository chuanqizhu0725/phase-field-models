#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

class Fields
{
public:
    int phaseNum;
    int Nx, Ny, Nz;
    Fields(int Nx0, int Ny0, int Nz0, int phaseNum0)
    {
        Nx = Nx0;
        Ny = Ny0;
        Nz = Nz0;
        phaseNum = phaseNum0;
    }

    double PhaseField()
    {
        // double pf = 0.1;
        // double pf[2];
        // pf[0] = 0.1;
        // pf[1] = 0.2;
        // return *pf;
        // cout << sizeof(pf) << endl;
        double pf2d[phaseNum][Nx][Ny];
        double pf3d[phaseNum][Nx][Ny][Nz];
        if (Nz == 0)
        {
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < phaseNum; k++)
                    {
                        pf2d[k][i][j] = 0.0;
                    }
                }
            }
            return pf2d;
        }
        else if (Nz > 0)
        {
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int l = 0; l < Nz; l++)
                    {
                        for (int k = 0; k < phaseNum; k++)
                        {
                            pf3d[k][i][j][l] = 0.0;
                        }
                    }
                }
            }
            return pf3d;
        }
    }
};