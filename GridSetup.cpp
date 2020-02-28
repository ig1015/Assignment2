#include <math.h>

void GridSetup (int * coord, int size, int Px, int Py, int Nx_g, int Ny_g, int& Nx, int& Ny)
{
    int Nx_sub = ceil( (1.0 * Nx_g / Px ) ) +2;  // 2 added for boundary points  
    int Ny_sub = ceil( (1.0*Ny_g / Py) ) + 2;
    int Nx_last = Nx_g - (Px-1) * (Nx_sub-2) +1;      // last domain takes over leftover points
    int Ny_last = Ny_g - (Py-1) * (Ny_sub-2) +1;
    
    if ( coord[0] == (Py -1)  )
        {
            Ny = Ny_last;
        }
        else
        {
            Ny = Ny_sub;
        }
        
        if ( coord[1] == (Px -1) ) 
        {
            Nx = Nx_last;
        }
        else
        {
            Nx = Nx_sub;
        }
        
        if (coord[0] == 0)
        {
            Ny -= 1;
        }
        if (coord[1] == 0)
        {
            Nx -= 1;
        }
}