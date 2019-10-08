#ifndef PARTTYPE_H
#define PARTTYPE_H

#include "pointers.h"
#include <string>
#include <vector>

/*
#include <stdlib.h>

#include "mpi.h"

#include <stdio.h>
#include <iostream>
#include <fstream>


#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
*/


#define MASTER 0

namespace MASKE_NS {
    
	class Particles : protected Pointers {
	public:
		
		 Particles(class MASKE *);
		~Particles();
		
        std::vector<std::string> typenames;
        int Nt, Np;
        //LAMMPS_NS::LAMMPS *lmp;
        
        void addtype(std::string);
        void printall();
		
	private:

        
	};
	
}

#endif
