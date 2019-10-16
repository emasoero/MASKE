#ifndef BLOCK_H
#define BLOCK_H

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
	
	class Block : protected Pointers {
	public:
        Block(class MASKE *);
		~Block();
    
        int me;     // id of the current processor (rank)
        
        int Nreg, Nlat; //number of recorded region and lattice blocks

        std::vector<std::string> LatNames;    // names of all lattice blocks input by user
        std::vector<std::string> RegNames;  // names of all region blocks input by user
        
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
        void printall();
        
        
	private:

	};
	
}

#endif
