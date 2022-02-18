#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

#define Nx 128

int main()
{
    double dx = 1.0;
    double dtime = 0.2;
    int nstep = 600;
    int nprint = 200;

    double temp[Nx];

    for (int i = 0; i < Nx; i++)
    {
        if (i > 44 && i < 84)
        {
            temp[i] = 1.0;
        }
        else
        {
            temp[i] = 0.0;
        }
    }

    for (int istep = 0; istep < nstep; istep++)
    {
        for (int i = 0; i < Nx; i++)
        {
            temp[i] = temp[i] + dtime * (temp[i + 1] + temp[i - 1] - 2.0 * temp[i]) / (dx * dx);
        }
        if (istep % nprint == 0)
        {
            ostringstream oss;
            oss << "data/1d" << istep << ".csv";
            ofstream MyFile(oss.str());
            for (int j = 0; j < Nx; j++)
            {
                MyFile << temp[j] << "\n";
            }
            // Close the file
            MyFile.close();
        }
    }
    return 0;
}