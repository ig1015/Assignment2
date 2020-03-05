#include "PoissonSolver.h"
#include "cstring"
#include "iostream"

#define F77NAME(x) x##_

extern "C" {
void F77NAME(dpbtrf) (
const char& UPLO,
const int& n, const int& kd,
const double* A, const int& ldab,
const int& info);
}

extern "C" {
void F77NAME(dpbtrs) (
const char& UPLO,
const int& n, const int& kd, const int& NRHS,
const double* A, const int& ldab,
double* b, const int& LDB,
const int& info);
}



PoissonSolver::PoissonSolver(int Nx_in, int Ny_in, double dx_in, double dy_in)
{
    Nx = Nx_in;
    Ny = Ny_in;
    dx = dx_in;
    dy = dy_in;
    Ny_r = Ny + 1;
    A = new double [(Nx)*(Ny)*(Ny+1)];
    memset(A, 0, (Nx-2)*(Ny-2)*(Ny-1)*sizeof(double));
}

PoissonSolver::~PoissonSolver()
{
}

void PoissonSolver::Matrix()
{
    for (int i = 0 ; i < (Nx)*(Ny) ; i++)
    {
        for (int j = 0 ; j < (Ny_r) ; j++)
        {
            if (j == (Ny_r-1) )
            {
                A[i*(Ny_r) + j] = 2 / dx /dx + 2/ dy /dy;
            }
            if ( (j ==(Ny_r -2)) && ( (i < 1) || i % (Ny) != 0 ) && i>0 )
            {
                A[i*(Ny_r) + j] = -1/ dy / dy;
            }
            if (j == 0 && i > (Ny-1))
            {
                A[i*(Ny_r) + j] = -1 / dx /dx;
            }
        }
    }
    
    /*for (int i = 0 ; i < (Nx-2)*(Ny-2) ; i++)
    {
        for (int j = 0 ; j < (Ny_r) ; j++)
        {
            if (j == (Ny_r-1) )
            {
                A[i*(Ny_r) + j] = 2 / dx /dx + 2/ dy /dy;
            }
            if ( (j ==(Ny_r -2)) && ( (i < 1) || i % (Ny-2) != 0 ) && i>0 )
            {
                A[i*(Ny_r) + j] = -1/ dy / dy;
            }
            if (j == 0 && i > (Ny-3))
            {
                A[i*(Ny_r) + j] = -1 / dx /dx;
            }
        }
    }*/
}

void PoissonSolver::Initialise()
{
    Matrix();
    Factorise();
}

void PoissonSolver::Integrate(double* RHS)
{
    /*cout << "A- ugabuga" << endl;
    for (int j = 0 ; j<(Ny_r) ; j++)
    {
        for (int i = 0 ; i<(Nx)*(Ny) ; i++)
        {
            cout.width(4);
            cout << A[i*Ny_r + j] << " ";
        }
        cout << endl;
    }*/
    
    // Solve with known factorisation
    info = 0;
    F77NAME(dpbtrs)('U', Nx*Ny, Ny, 1, A, Ny+1, RHS, Nx*Ny, info);
    
    
}

void PoissonSolver::Factorise()
{
    info = 0;
    F77NAME(dpbtrf)('U', Nx*Ny, Ny, A, Ny+1, info);
}