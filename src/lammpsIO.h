#ifndef LAMMPSIO_H
#define LAMMPSIO_H

#include "pointers.h"
#include "mpi.h"
#include <string>
#include <sstream>
/*
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"



#define MASTER 0

namespace MASKE_NS {
    
	class LammpsIO : protected Pointers {
	public:
		
		 LammpsIO(class MASKE *);
		~LammpsIO();
		
        LAMMPS_NS::LAMMPS *lmp;
        
        void create();  //Create a new Lammps instance
        void printall();
        void lammpsdo(std::string); //passes a command to lammps
        
        std::string units;
        std::string atomstyle;
        std::string lmpThSt;    // string recording the themo style from input: used to add and evaluate tem computes when needed in source code of lammps
        double timestep;
        std::stringstream ss;
        bool lammps_active; // makes the current processor if lammps has been already activated for its subcomm
        //std::string tdump_fname;  // temporary dump file name, used by maske to evaluate computes etc for output
        
	private:
		//std::string fname;              // open inputcprs file
        int me;     // id of the current processor (rank)
		
	};
	
}

#endif
