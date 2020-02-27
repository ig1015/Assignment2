#include "LidDrivenCavity.h"
#include <cstring>
#include <iostream>

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx = nx;
    Ny = ny;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re = re;
}

void LidDrivenCavity::Initialise()
{
    v = new double [Nx*Ny];
    memset(v,0,Nx*Ny*sizeof(double));
    s = new double [Nx*Ny];
    memset(s,0,Nx*Ny*sizeof(double));
}

void LidDrivenCavity::Integrate()
{
}

double LidDrivenCavity::get_dt()
{
    return dt;
}

double LidDrivenCavity::get_Nx()
{
    return Nx;
}

double LidDrivenCavity::get_Ny()
{
    return Ny;
}

void LidDrivenCavity::printv()
{
    for (int i = 0 ; i<Nx ; i++)
    {
        for (int j = 0 ; j<Ny ; j++)
        {
            cout << v[i*Ny + j] << " ";
        }
        cout << endl;
    }
}

void LidDrivenCavity::SetCartesianCoordinates(int y , int x)
{
    y_c = y;
    x_c = x;
}