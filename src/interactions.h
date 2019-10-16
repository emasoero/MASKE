#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include "pointers.h"
#include <string>
#include <vector>
//#include <stdlib.h>
/*#include "mpi.h"

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
    
	class Interactions : protected Pointers {
	public:
		
        Interactions(class MASKE *);
		~Interactions();
		
        std::string stylestr;
        std::vector<std::string> Mcoeff;
        int Nc;
        
        //LAMMPS_NS::LAMMPS *lmp;
        
        void printall();
        void writestyle(std::string);
        void addpcoeff(std::string,std::string,std::string);
		
	private:

        
	};
	
}

#endif
