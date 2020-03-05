#pragma once

#include <string>
using namespace std;

class PoissonSolver
{
public:
    PoissonSolver(int Nx_in, int Ny_in, double dx_in, double dy_in);
    ~PoissonSolver();

    void Matrix();
    void Factorise();

    void Initialise();
    void Integrate(double* RHS);

    // Add any other public functions

private:
    double* b = nullptr;
    double* A = nullptr;
    
    int Nx;
    int Ny;
    int Ny_r;
    double dx;
    double dy;
    int info;
    
    
    
    
};
