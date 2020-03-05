#include "LidDrivenCavity.h"
#include <cstring>
#include <iostream>

#define F77NAME(x) x##_

extern "C" {
void F77NAME(dpbsv) (
const char& UPLO,
const int& n, const int& kd,
const int& nrhs,
const double* A, const int& ldab,
const double* y, const int& ldb,
const int& info);
}

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
    vNew = new double [(Nx-2)*(Ny-2)];
    memset(vNew,0,(Nx-2)*(Ny-2)*sizeof(double));
    s = new double [Nx*Ny];
    memset(s,0,Nx*Ny*sizeof(double));
    A = new double [(Nx-2)*(Ny-2)*(Ny-1)];
    memset(A, 0, (Nx-2)*(Ny-2)*(Ny-1)*sizeof(double));
    fill_n(s , Nx*Ny , rank);               // using initialise to rank to check sending
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
    int k, info, Ny_r;
    double timer = 0;
    while (timer < T)
    {
        // Obtain vorticity at boundaries at time t
    Boundary();
    
    // why is this here 
    //mySend(s);
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Calculate vorticity at time t
    for (int i=1 ; i < Nx-1 ; i++ )         //i index restricted from second to penultimate node 
    {
        for (int j=1 ; j<Ny-1 ; j++)
        {
            v[i*Ny + j] = - ( (s[(i+1)*Ny+j] - 2.0 * s[i*Ny + j] + s[(i-1)*Ny+j] ) / dx /dx 
                + ( s[i*Ny + j+1] - 2.0 * s[i*Ny + j] + s[i*Ny + j-1] ) / dy / dy );
        }
    }
    

    
    mySend(v);
    
    //////////////////////////////////////////////////////////////////////////
    k = 0;
    for (int i=1 ; i < Nx-1 ; i++ )         //i index restricted from second to penultimate node 
    {
        for (int j=1 ; j<Ny-1 ; j++)
        {
            vNew[k] = v[i*Ny+j] + dt * (s[i*Ny+j+1] - s[i*Ny+j-1])  
                * (v[(i+1)*Ny+j] - v[(i-1)*Ny+j]) / 4.0 / dx /dy  
                - dt * (s[(i+1)*Ny+j] - s[(i-1)*Ny+j]) * (v[i*Ny+j+1] 
                - v[i*Ny+j-1]) / 4.0 / dx /dy 
                + dt / Re * ( (v[(i+1)*Ny+j] - 2.0 * v[i*Ny + j] + v[(i-1)*Ny+j] ) / dx / dx 
                + ( v[i*Ny + j+1] - 2.0 * v[i*Ny + j] + v[i*Ny + j-1] ) / dy / dy );
                
                k++;
                
                
   /*vNew[(i-1)*(Ny-2) + (j-1)] = v[i*Ny + j] + dt * (  1.0 / Re * ( (v[(i+1)*Ny+j] - 2 * v[i*Ny + j] + v[(i-1)*Ny+j] ) / (dx * dx) 
                + ( v[i*Ny + j+1] - 2 * v[i*Ny + j] + v[i*Ny + j-1] ) / (dy *dy)  ) 
                    - ( ( ( s[i*Ny + j-1] - s[i*Ny + j+1] ) / 2 / dy ) * ( ( v[(i+1)*Ny+j] - v[(i-1)*Ny+j] ) / 2 /dx ) )
                        + ( ( ( v[i*Ny + j-1] - v[i*Ny + j+1] ) / 2 / dy ) * ( ( s[(i+1)*Ny+j] - s[(i-1)*Ny+j] ) / 2 /dx ) ) );*/
        }
    }
    
    // set new vorticity value
    
    k = 0;
    for (int i=1 ; i < Nx-1 ; i++ )         //i index restricted from second to penultimate node 
    {
        for (int j=1 ; j<Ny-1 ; j++)
        {
            v[i*Ny + j] = vNew[k];
            k++;
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    // build coefficient matrix
    
    Ny_r = Ny-1;
    
    memset(A, 0, (Nx-2)*(Ny-2)*(Ny-1)*sizeof(double));
    
    for (int i = 0 ; i < (Nx-2)*(Ny-2) ; i++)
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
    }
    
    cout << "A" << endl;
    for (int j = 0 ; j<(Ny_r) ; j++)
    {
        for (int i = 0 ; i<(Nx-2)*(Ny-2) ; i++)
        {
            cout.width(4);
            cout << A[i*Ny_r + j] << " ";
        }
        cout << endl;
    }
    
    

    ///////////////////////////////////////////////////////////////////////////////////////
    
    // Solve equation
    
    info = 0;

    F77NAME(dpbsv)('U', (Ny-2) * (Nx-2), (Ny-2), 1, A, (Ny-1), vNew, (Ny-2)*(Nx-2), info);
    
    k = 0;
    for (int i=1 ; i < Nx-1 ; i++ )         //i index restricted from second to penultimate node 
    {
        for (int j=1 ; j<Ny-1 ; j++)
        {
            s[i*Ny + j] = vNew[k];
            k++;
        }
    }
    
    // Check if parameters are correct
    
    //cout << dx << endl;
    //cout << dy << endl;
    
    
    timer += dt;
    }
    
    
    
}

void LidDrivenCavity::mySend(double * ar)
{
     if(neigh[0] != -2)
    {
        // copy 2nd entry of each column (i.e. 2nd row) 
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
            vert_comm_t[i] = ar[Ny*(i+1)+1];
        }
    }
    
    if(neigh[1] != -2)
    {
        // copy penultimate column
        for (int i=0 ; i < (Ny-2) ; i++)
        {
            hori_comm_r[i] = ar[(Nx-2)*Ny+i+1];
        }
    }
    
    // If there exists a neighbour north
    if(neigh[2] != -2)
    {
        // copy penultimate entry of each column (i.e. 2nd row) 
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
            vert_comm_b[i] = ar[Ny*(i+1)+(Ny-2)];
        }
    }
    
    if (neigh[3] != -2)
    {
        // copy 2nd column
        for (int i = 0 ; i < (Ny -2) ; i++)
        {
            hori_comm_l[i] = ar[Ny + i + 1];
        }
    }
   // Send out all the relevant boundaries
    //send to top neighbour
    if(neigh[0] != -2)
    {
        MPI_Send(vert_comm_t, Nx-2, MPI_DOUBLE, neigh[0], 0, MPI_COMM_WORLD);
    }
    // send to neighbout on the right
    if(neigh[1] != -2)
    {
        MPI_Send(hori_comm_r, Ny-2, MPI_DOUBLE, neigh[1], 0, MPI_COMM_WORLD);
    }
    // send to bottom neighbour
    if(neigh[2] != -2)
    {
        MPI_Send(vert_comm_t, Nx-2, MPI_DOUBLE, neigh[2], 0, MPI_COMM_WORLD);
    }
    // send to left neighbour
    if (neigh[3] != -2)
    {
        MPI_Send(hori_comm_l, Ny-2, MPI_DOUBLE, neigh[3], 0, MPI_COMM_WORLD);
    }
    
    // Receive messages
    // receive from top neighbour
    if (neigh[0] != -2)
    {
        MPI_Recv(vert_comm_t, Nx-2, MPI_DOUBLE, neigh[0], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
             ar[Ny*(i+1)] = vert_comm_t[i];
        }
    }
    //Receive message from right neighbour
    if (neigh[1] != -2)
    {
        MPI_Recv(hori_comm_r, Ny-2, MPI_DOUBLE, neigh[1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i=0 ; i < (Ny-2) ; i++)
        {
             ar[(Nx-1)*Ny+i+1] = hori_comm_r[i];
        }
    }
    //receive message from bottom neighbour
    if (neigh[2] != -2)
    {
        MPI_Recv(vert_comm_b, Nx-2, MPI_DOUBLE, neigh[2], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
             ar[Ny*(i+1)+(Ny-1)] = vert_comm_b[i];
        }
    }
    //Receive message from left neighbour
    if (neigh[3] != -2)
    {
        MPI_Recv(hori_comm_l, Ny-2, MPI_DOUBLE, neigh[3], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < (Ny -2) ; i++)
        {
            ar[i + 1] = hori_comm_l[i];
        }
    }
}


void LidDrivenCavity::Boundary()
{
    
    
    /*// top boundary
    if (neigh[0] == -2)
    {
        for (int i = 0; i<Nx; i++)
        {
            // assume U=1, I*Ny -> first element in each column (Data is saved column wise)
            v[i*Ny + (Ny - 1)] =   2.0 / dy/dy * (s[i*Ny + (Ny - 1)] - s[i*Ny + (Ny - 2)]) - 2.0 / dy;
        }
    }
    
    // bottom boundary
    if (neigh[2] == -2)
    {
        for (int i = 0; i<Nx; i++)
        {
            v[i*Ny] = 2.0/ (dy*dy) * ( s[i*Ny] - s[i*Ny + 1] ) ;
        }
    }
    
    // right boundary
    if (neigh[1] == -2)
    {
        for (int i = 0 ; i < Ny ; i++)
        {
            v[(Nx-1)*Ny+i] = 2.0 / dx / dx * (s[(Nx-1)*Ny+i] - s[(Nx-2)*Ny+i]) ;
        }
    }
    
    if (neigh[3] == -2)
    {
        for (int i = 0 ; i < Ny ; i++)
        {
            v[i] = 2.0 / dx / dx * (s[i] - s[Ny+i]) ;
        }
    }*/
     
     
    
    
      // top boundary
    if (neigh[0] == -2)
    {
        for (int i = 0; i<Nx; i++)
        {
            // assume U=1, I*Ny -> first element in each column (Data is saved column wise)
            v[i*Ny] = 2.0/ (dy*dy) * ( s[i*Ny] - s[i*Ny + 1] ) - 2.0 / dy ;
        }
    }
    
    // right boundary
    if (neigh[1] == -2)
    {
        for (int i = 0 ; i < Ny ; i++)
        {
            v[(Nx-1)*Ny+i] = 2.0 / (dx * dx) * (s[(Nx-1)*Ny+i] - s[(Nx-2)*Ny+i]) ;
        }
    }
    
    // bottom boundary
    if (neigh[2] == -2)
    {
        for (int i = 0; i<Nx; i++)
        {
            v[i*Ny + (Ny - 1)] =   2.0 / (dy*dy) * (s[i*Ny + (Ny - 1)] - s[i*Ny + (Ny - 2)]) ;
        }
    }
    
    // left boundary
    if (neigh[3] == -2)
    {
        for (int i = 0 ; i < Ny ; i++)
        {
            v[i] = 2 / dx / dx * (s[i] - s[Ny+i]) ;
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
            cout.width(10);
            cout << v[i + j*Ny] << " ";
        }
        cout << endl;
    }
}

void LidDrivenCavity::prints()
{
    for (int i = 0 ; i<Ny ; i++)
    {
        for (int j = 0 ; j<Nx ; j++)
        {
            cout.width(10);
            cout << s[i + j*Ny] << " ";
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


// All the MPI stuff now in function

////////////////////////////////////////////////////////////////////////////
    // CHECK THE INDEXING, POTENTIAL SOURC OF ERROR 
    ////////////////////////////////////////////////////////////////////////////
    
    // If there exists a neighbour north
    /*if(neigh[0] != -2)
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
            vert_comm_b[i] = s[Ny*(i+1)+(Ny-2)];
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
    
    
    ////////////////////////////////////////////////////////////////////////
    
    // Send out all the relevant boundaries
    //send to top neighbour
    if(neigh[0] != -2)
    {
        MPI_Send(vert_comm_t, Nx-2, MPI_DOUBLE, neigh[0], 0, MPI_COMM_WORLD);
    }
    // send to neighbout on the right
    if(neigh[1] != -2)
    {
        MPI_Send(hori_comm_r, Ny-2, MPI_DOUBLE, neigh[1], 0, MPI_COMM_WORLD);
    }
    // send to bottom neighbour
    if(neigh[2] != -2)
    {
        MPI_Send(vert_comm_t, Nx-2, MPI_DOUBLE, neigh[2], 0, MPI_COMM_WORLD);
    }
    // send to left neighbour
    if (neigh[3] != -2)
    {
        MPI_Send(hori_comm_l, Ny-2, MPI_DOUBLE, neigh[3], 0, MPI_COMM_WORLD);
    }
    
    // Receive messages
    // receive from top neighbour
    if (neigh[0] != -2)
    {
        MPI_Recv(vert_comm_t, Nx-2, MPI_DOUBLE, neigh[0], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
             s[Ny*(i+1)] = vert_comm_t[i];
        }
    }
    //Receive message from right neighbour
    if (neigh[1] != -2)
    {
        MPI_Recv(hori_comm_r, Ny-2, MPI_DOUBLE, neigh[1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i=0 ; i < (Ny-2) ; i++)
        {
             s[(Nx-1)*Ny+i+1] = hori_comm_r[i];
        }
    }
    //receive message from bottom neighbour
    if (neigh[2] != -2)
    {
        MPI_Recv(vert_comm_b, Nx-2, MPI_DOUBLE, neigh[2], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < (Nx-2) ; i++)
        {
             s[Ny*(i+1)+(Ny-1)] = vert_comm_b[i];
        }
    }
    //Receive message from left neighbour
    if (neigh[3] != -2)
    {
        MPI_Recv(hori_comm_l, Ny-2, MPI_DOUBLE, neigh[3], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < (Ny -2) ; i++)
        {
            s[i + 1] = hori_comm_l[i];
        }
    }
    */
    // Merge changes
    
    
    