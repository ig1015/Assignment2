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
    vNew = new double [Nx*Ny];
    memset(vNew,0,Nx*Ny*sizeof(double));
    s = new double [Nx*Ny];
    memset(s,0,Nx*Ny*sizeof(double));
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    vert_comm_t = new double [Nx-2];             // the corner nodes are not required for communication
    memset(vert_comm_t, 0, (Nx-2)*sizeof(double));
    vert_comm_b = new double [Nx-2];
    memset(vert_comm_b, 0, (Nx-2)*sizeof(double));
    hori_comm_r = new double [Ny -2];
    memset(hori_comm_r, 0, (Ny-2)*sizeof(double));
    hori_comm_l = new double [Ny -2];
    memset(hori_comm_l, 0, (Ny-2)*sizeof(double));
}

void LidDrivenCavity::Integrate()
{
    // Obtain vorticity at boundaries at time t
    Boundary();
    
    ////////////////////////////////////////////////////////////////////////////
    // CHECK THE INDEXING, POTENTIAL SOURC OF ERROR 
    ////////////////////////////////////////////////////////////////////////////
    
    // If there exists a neighbour north
    if(neigh[0] != -2)
    {
        // copy 2nd entry of each column (i.e. 2nd row) 
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
            vert_comm_t[i] = s[Ny*(i+1)+1];
        }
    }
    
    if(neigh[1] != -2)
    {
        // copy penultimate column
        for (int i=0 ; i < (Ny-2) ; i++)
        {
            hori_comm_r[i] = s[(Nx-2)*Ny+i+1];
        }
    }
    
    // If there exists a neighbour north
    if(neigh[2] != -2)
    {
        // copy penultimate entry of each column (i.e. 2nd row) 
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
            vert_comm_t[i] = s[Ny*(i+1)+(Ny-2)];
        }
    }
    
    if (neigh[3] != -2)
    {
        // copy 2nd column
        for (int i = 0 ; i < (Ny -2) ; i++)
        {
            hori_comm_l[i] = s[Ny + i + 1];
        }
    }
    
    // Calculate vorticity at time t
    for (int i=1 ; i < Nx-1 ; i++ )         //i index restricted from second to penultimate node 
    {
        for (int j=1 ; j<Ny ; j++)
        {
            v[i*Ny + j] = (s[(i+1)*Ny+j] - 2 * s[i*Ny + j] + s[(i-1)*Ny+j] ) / (dx * dx) 
                + ( s[i*Ny + j+1] - 2 * s[i*Ny + j] + s[i*Ny + j+1] ) / (dy *dy);
        }
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    // Send out all the relevant boundaries
    if(neigh[0] != -2)
    {
    }
    
    if(neigh[1] != -2)
    {
    }
    
    // If there exists a neighbour north
    if(neigh[2] != -2)
    {
    }
    
    if (neigh[3] != -2)
    {
    }
    
}

void LidDrivenCavity::Boundary()
{
    // top boundary
    if (neigh[0] == -2)
    {
        for (int i = 0; i<Nx; i++)
        {
            // assume U=1, I*Ny -> first element in each column (Data is saved column wise)
            v[i*Ny] = 2/ (dy*dy) * ( s[i*Ny] - s[i*Ny + 1] ) - 2 / dy;
        }
    }
    
    // right boundary
    if (neigh[1] == -2)
    {
        for (int i = 0 ; i < Ny ; i++)
        {
            v[(Nx-1)*Ny+i] = 2 / (dx * dx) * (s[(Nx-1)*Ny+i] - s[(Nx-2)*Ny+i]);
        }
    }
    
    // bottom boundary
    if (neigh[2] == -2)
    {
        for (int i = 0; i<Nx; i++)
        {
            v[i*Ny + (Ny - 1)] = 2 / (dy*dy) * (s[i*Ny + (Ny - 1)] - s[i*Ny + (Ny - 2)]);
        }
    }
    
    // left boundary
    if (neigh[3] == -2)
    {
        for (int i = 0 ; i < Ny ; i++)
        {
            v[i] = 2 / (dx * dx) * (s[i] - s[Ny+i]);
        }
    }
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
    for (int i = 0 ; i<Ny ; i++)
    {
        for (int j = 0 ; j<Nx ; j++)
        {
            cout << v[i + j*Ny] << " ";
        }
        cout << endl;
    }
}

void LidDrivenCavity::SetCartesianCoordinates(int y , int x)
{
    y_c = y;
    x_c = x;
}

void LidDrivenCavity::SetNeighbours(int *p, int size)
{
    // Convention is such that neighbours are defined: North, East, South, West
    for (int i=0 ; i<size ; i++ )
    {
            neigh[i] = p[i];
    }
}

void LidDrivenCavity::getNeighbours(int*p)
{
   for (int i=0 ; i<4 ; i++ )
    {
            p[i] = neigh[i];
    } 
}

void LidDrivenCavity::SetRank(int processRank)
{
    rank = processRank;
}

void LidDrivenCavity::SetComms(MPI_Comm group)
{
    mycomms = group;
}