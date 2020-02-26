#include <iostream>
#include <boost/program_options.hpp>
#include <math.h>
using namespace std;

namespace po = boost::program_options;

#include "LidDrivenCavity.h"
#include <mpi.h> 

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
        

            // define subdomain grid sizes
            
            int Nx_sub = ceil( (1.0 * vm["Nx"].as<int>()/vm["Px"].as<int>()) );       
            cout << Nx_sub << endl;
            
            int Ny_sub = ceil( (1.0*vm["Ny"].as<int>()/vm["Py"].as<int>()) );
            
            int Nx_last = vm["Nx"].as<int>() - (vm["Px"].as<int>()-1) * Nx_sub;      // last domain takes over leftover points
            cout << "last domain takes " << Nx_last << endl;
            
            int Ny_last = vm["Ny"].as<int>() - (vm["Py"].as<int>()-1) * Ny_sub;
        
    
    

    MPI_Finalize();
    
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