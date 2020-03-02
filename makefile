default: myprog

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp LidDrivenCavity.h GridSetup.h
	mpicxx -std=c++11 -Wall -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp 
	
LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h
	mpicxx -std=c++11 -Wall -o LidDrivenCavity.o -c LidDrivenCavity.cpp
	
GridSetup.o: GridSetup.cpp GridSetup.h
	mpicxx -std=c++11 -Wall -o GridSetup.o -c GridSetup.cpp
	
myprog: LidDrivenCavitySolver.o LidDrivenCavity.o GridSetup.o
	mpicxx -o myprog LidDrivenCavitySolver.o LidDrivenCavity.o GridSetup.o -lboost_program_options
	
.PHONY: clean

clean:
	rm -f *.o myprog