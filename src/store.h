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

        /*
        std::vector<std::string> tempSvec;  // temp vector to store all lammps commands associated to block before pushing them back into vector of vectors
        
        
        std::vector<std::vector<std::string> > LatLammps; // collaction of all lammps command lines in each lattice block
        std::vector<std::vector<std::string> > RegLammps; // collaction of all lammps command lines in each region block
        
        bool MinModFound, MinCmdFound; // flags indicating whether Minimiser and Min Modifyer are found (as they should) in lattice block
        bool LatFound; // flag indicating whether "lattice" command is found (as it should) in lattice block
        bool RegFound; // flag indicating whether "region" command is found (as it should) in region block

        bool BlockFound;    // flag indicating whether another block command is found within a block. This should not happen: if it happens it means that a previous block was not correclty closed via the command end_block
        
        std::vector<std::string> LatMinMod;     // minimization style of each lattice
        std::vector<std::string> LatMinCmd;    // minimization command of each lattice
        
        
        
        void add_lattice(std::string);
        void add_region(std::string);
        void add_Latline(std::string);
        void add_Regline(std::string);
        */
        void add(std::string);
        void printall();
        
        
	private:

	};
	
}

#endif
