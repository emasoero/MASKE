#ifndef SIMBOX_H
#define SIMBOX_H

#include "pointers.h"
#include <string>
#include <stdlib.h>
#include <vector>
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
    
	class Simbox : protected Pointers {
	public:
		
		 Simbox(class MASKE *);
		~Simbox();
		
        std::vector<std::string> boundary;
        
        double xyzlo[3], xyzhi[3], xyztri[3];
        bool triclinic;
        
        
        double** mtri; //  matrix to transform from normalized 01-01-01 to triclinic LAMMPS box, according to the equation   Vtri =  mtri  *  V01  + xyzlo, where  Vtri is a column vector in the triclinic box, V01 is the corresponding column vector in the normalized box, * indicates the row-by-colum product, and xyzlo is the origin of the triclinic box according to LAMMPS
        
        double** mtriINV; // matrix to transform from triclinic LAMMPS box to normalized 01-01-01 box, according to the equation  V01 =  mtriINV  *  (Vtri - xyzlo)
        
        //LAMMPS_NS::LAMMPS *lmp;
        
        void printall();
        void computeTM(); // compute transformation matrices to go from triclinic to normalised space and backs
        void MVprod(double**, double*);
        void N2Ttrans(double**, double*, double*);
        void T2Ntrans(double**, double*, double*);
        void MatInv(double**);
		
	private:

        
	};
	
}

#endif
