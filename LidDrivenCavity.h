#pragma once

#include <string>
#include <mpi.h> 

using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetCartesianCoordinates(int y , int x);
    void SetNeighbours(int * p, int size);
    void SetRank(int rank);
    // Will this work, it doesn't throw an error but how do you pass it to the object
    void SetComms(MPI_Comm group);

    void Initialise();
    void Integrate();
    void Boundary();
    
    double get_dt();
    double get_Nx();
    double get_Ny();
    void printv();
    void prints();
    void getNeighbours(int*p);

    // Add any other public functions

private:
    double* v = nullptr;
    double* vNew = nullptr;
    double* s = nullptr;
    int neigh[4];
    double* vert_comm_t = nullptr;
    double* hori_comm_r = nullptr;
    double* vert_comm_b = nullptr;
    double* hori_comm_l = nullptr;

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    int x_c;
    int y_c;
    int rank;
    double dx;
    double dy;
    MPI_Comm mycomms;
};
