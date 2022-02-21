#include <iostream>
using namespace std;

class Domain
{
public:
    double Lx, Ly, Lz;
    int Nx, Ny, Nz;
    double dx, dy, dz;
    Domain(double Lx0, double Ly0, double Lz0, int Nx0, int Ny0, int Nz0, double dx0, double dy0, double dz0)
    {
        Lx = Lx0;
        Ly = Ly0;
        Lz = Lz0;
        Nx = Nx0;
        Ny = Ny0;
        Nz = Nz0;
        dx = dx0;
        dy = dy0;
        dz = dz0;
    }

    void about()
    {
        cout << Lx << "\n";
        cout << Ly << "\n";
        cout << Lz << "\n";
        cout << Nx << "\n";
        cout << Ny << "\n";
        cout << Nz << "\n";
        cout << dx << "\n";
        cout << dy << "\n";
        cout << dz << "\n";
    };
};
