#ifndef STORE_H
#define STORE_H

#include "pointers.h"
#include <string>
#include <vector>

//#include "mpi.h"
//#include <string>
//#include <vector>
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#define MASTER 0

namespace MASKE_NS {
	
	class Store : protected Pointers {
	public:
        Store(class MASKE *);
		~Store();
    
        int me;     // id of the current processor (rank)
        int Nreg, Nlat, Nmin;   // number of stored regions, lattices, and minimisers
        
        std::vector<std::string> RegNames;  // names of all regions stored by user
        std::vector<std::string> LatNames;    // names of all lattices stored by user
        std::vector<std::string> MinNames;  // names of all minimisers stored by user
	std::vector<std::string> MulNames;

        std::vector<std::string> RegCmd;  // vector with all region commants stored by user
        std::vector<std::string> LatCmd;    // same for lattices
        std::vector<std::string> LatDVnames;   // see below
        std::vector<std::string> LatDVcmd; // vectors with LAMMPS variable names and commands to compute DV associated to a lattice point, to be used in fix_nucleate sampling
        std::vector<std::string> MinCmd;  // same for minimize
        std::vector<std::string> MinTstep;  // time step for minimizer
        std::vector<std::string> MinModCmd;  // same for min_modify

        std::vector<std::vector<std::string>> MulCmd;

        void add(std::string);
        void printall();
        
        
	private:

	};
	
}

#endif
