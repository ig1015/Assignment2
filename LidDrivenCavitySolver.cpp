#include <iostream>
#include <boost/program_options.hpp>
#include <math.h>
using namespace std;

namespace po = boost::program_options;

#include "LidDrivenCavity.h"
#include <mpi.h> 
#include "GridSetup.h"

int main(int ac, char **av)
{
    // Getting user inputs using boost
     
    MPI_Init(&ac, &av);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
     
        
        // Try block used to ensure user provides all the required inputs.
        // Within the try block the values for class data members are assigned. Checks are carried out to ensure the class
        // variables are compatible. Once all data members are assigned the vecotr map created goes out of scope. 
        po::variables_map vm;
        try 
        {
            po::options_description desc("Allowed options");
            desc.add_options()
                ("help", "produce help message")
                ("dt", po::value<double>(), "set timestep")
                ("T", po::value<double>(), "set total simulation time")
                ("Nx", po::value<int>(), "set number of points in x-direction")
                ("Ny", po::value<int>(), "set number of points in y-direction")
                ("Re", po::value<double>(), "set Reynold's number")
                ("Lx", po::value<double>(), "set x-length of domain")
                ("Ly", po::value<double>(), "set Reynold's number")
                ("Px", po::value<int>(), "set number of x partitions")
                ("Py", po::value<int>(), "set number of y partitions")
            ;

            po::store(po::parse_command_line(ac, av, desc), vm);
            po::notify(vm);  
            
            if (rank == 0)
            {
                if (vm.count("help")) {
                    cout << desc << "\n";
                    return 0;
                }

                if (vm.count("dt")) {
                    cout << "dt was set to " 
                         << vm["dt"].as<double>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("dt not provided");
                }
                if (vm.count("T")) {
                    cout << "T was set to " 
                         << vm["T"].as<double>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("T not provided");
                }
                if (vm.count("Nx")) {
                    cout << "Nx was set to " 
                         << vm["Nx"].as<int>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("Nx not provided");
                }
                if (vm.count("Ny")) {
                    cout << "Ny was set to " 
                         << vm["Ny"].as<int>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("Ny not provided");
                }
                if (vm.count("Re")) {
                    cout << "Re was set to " 
                         << vm["Re"].as<double>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("Re not provided");
                }
                if (vm.count("Lx")) {
                    cout << "Lx was set to " 
                         << vm["Lx"].as<double>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("Re not provided");
                }
                if (vm.count("Ly")) {
                    cout << "Ly was set to " 
                         << vm["Ly"].as<double>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("Re not provided");
                }
                if (vm.count("Px")) {
                    cout << "Px was set to " 
                         << vm["Px"].as<int>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("Px not provided");
                }
                if (vm.count("Py")) {
                    cout << "Py was set to " 
                         << vm["Py"].as<int>() << ".\n";
                }
                else
                {
                    throw std::out_of_range("Py not provided");
                }
                if (size != vm["Px"].as<int>() * vm["Py"].as<int>())
                {
                    throw std::out_of_range("Px * Py must equal np");
                }
            }
        }
        catch(exception& e) {
            cerr << "error: " << e.what() << "\n";
            MPI_Finalize();
            return 1;
        }
        catch(...) {
            cerr << "Exception of unknown type!\n";
        }
        
        if (rank == 0)
        {
                //LidDrivenCavity A; 
            LidDrivenCavity* solverMaster = new LidDrivenCavity();
                // Set class properties
            solverMaster -> SetTimeStep(vm["dt"].as<double>());
            solverMaster -> SetFinalTime(vm["T"].as<double>());
            solverMaster -> SetGridSize(vm["Nx"].as<int>(), vm["Ny"].as<int>());
            solverMaster -> SetReynoldsNumber(vm["Re"].as<double>());
            solverMaster -> SetDomainSize(vm["Lx"].as<double>(), vm["Ly"].as<double>());
           
            double test = solverMaster -> get_dt();
            cout << test << endl;
            
            solverMaster->Initialise();
            
            solverMaster -> printv();
        } 
        
        int Px = vm["Px"].as<int>();
        int Py = vm["Py"].as<int>();
        double dx = vm["Lx"].as<double>() / ( vm["Nx"].as<int>() - 1 );
        double dy = vm["Ly"].as<double>() / ( vm["Ny"].as<int>() - 1 );
        
        // Creat object for each subdomain - addd parameters which are common to all subdomains
        LidDrivenCavity* solver = new LidDrivenCavity();
        solver -> SetTimeStep(vm["dt"].as<double>());
        solver -> SetFinalTime(vm["T"].as<double>());
        solver -> SetReynoldsNumber(vm["Re"].as<double>());
        
       
        
        // Create cartesian grid
        MPI_Comm mygrid;
        const int dims = 2;
        int sizes[dims] = {vm["Py"].as<int>(), vm["Px"].as<int>()};
        int periods[dims] = {0,0};
        MPI_Cart_create(MPI_COMM_WORLD, dims, sizes , periods, 0, &mygrid);
    
        // get coordinates for each rank
        int coord[dims];
        MPI_Cart_coords(mygrid, rank, dims, coord);
        
        //get sizes of each domain
        /*int Nx_sub = ceil( (1.0 * vm["Nx"].as<int>()/vm["Px"].as<int>()) ) +2;  // 2 added for boundary points  
        int Ny_sub = ceil( (1.0*vm["Ny"].as<int>()/vm["Py"].as<int>()) ) + 2;
        int Nx_last = vm["Nx"].as<int>() - (vm["Px"].as<int>()-1) * (Nx_sub-2) +1;      // last domain takes over leftover points
        int Ny_last = vm["Ny"].as<int>() - (vm["Py"].as<int>()-1) * (Ny_sub-2) +1;*/
        
        
        
        // add elements depending on the position in the cartesian grid
        int Nx;
        int Ny;
        
        /*if ( coord[0] == (Py -1)  )
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
        }*/
        
    
        GridSetup(coord,3,Px,Py,vm["Nx"].as<int>(), vm["Ny"].as<int>(),  Nx,Ny);
        
        solver -> SetGridSize(Nx,Ny);
        solver -> SetDomainSize((Nx-1)*dx, (Ny-1)*dy);
        solver -> SetCartesianCoordinates(coord[0], coord[1]);
        
        MPI_Finalize();
        
        cout << "coordinates: " << coord[0] << ", " << coord[1] << "size" << Nx << ":" << Ny <<endl;
        
                

         /*// this ends up meaning that if we use a really course domain and loads of processes it will fail (could check and if that's the case 
            // use floor instead!!
            // define subdomain grid sizes
            
            
            int Nx_sub = ceil( (1.0 * vm["Nx"].as<int>()/vm["Px"].as<int>()) );       
            cout << Nx_sub << endl;
            
            int Ny_sub = ceil( (1.0*vm["Ny"].as<int>()/vm["Py"].as<int>()) );
            
            int Nx_last = vm["Nx"].as<int>() - (vm["Px"].as<int>()-1) * Nx_sub;      // last domain takes over leftover points
            cout << "last domain takes " << Nx_last << endl;
            
            int Ny_last = vm["Ny"].as<int>() - (vm["Py"].as<int>()-1) * Ny_sub;*/
            
            
            
            //double dx = vm["Lx"].as<double>() / ( vm["Nx"].as<int>() - 1 );
            //double dy = vm["Ly"].as<double>() / ( vm["Ny"].as<int>() - 1 );
            
            
       /*      // Asign appropriate sizes to subdomains
        if ( coords[0] != ( vm["Px"].as<int>() - 1 ) && ( coords[1] != vm["Py"].as<int>() - 1 ) )
            {
                solver -> SetGridSize(Nx_sub,Ny_sub);
                cout << "hi" << rank << endl; 
            }
        else if ( (rank >= vm["Py"].as<int>() * (vm["Px"].as<int>() - 1)) && ((rank+1) % vm["Py"].as<int>() != 0) )
            {
                solver -> SetGridSize(Nx_last,Ny_sub);
                cout << "dick" << rank << endl; 
            }
        else if ( (rank < vm["Py"].as<int>() * (vm["Px"].as<int>() - 1) ) && ((rank+1) % vm["Py"].as<int>() == 0) )
            {
                solver -> SetGridSize(Nx_sub,Ny_last);
                cout << "dog" << rank << endl;
            }
        else
            {
                solver -> SetGridSize(Nx_last,Ny_last);
                cout << "cat" << rank << endl;
            } */
            

            
            
            
        
    
    

    
    
    // Chris stuff
    /*// Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    // ...
    solver->Initialise();

    // Run the solver
    solver->Integrate();*/

	return 0;
}