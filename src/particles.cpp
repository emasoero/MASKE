#include "particles.h"
//#include "universe.h"
//#include "DTnucleate.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "clinker.h"
#include "pyramids.h"
#include "hydration.h"
#include "memory.h"
#include "csh.h"
#include "c2s.h"
#include "c3s.h"
#include "water.h"
#include <string.h>
*/

using namespace MASKE_NS;


Particles::Particles(MASKE *maske) : Pointers(maske)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    /*MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating lammpsIO class\n");

    lammps_active=false;
    units = "real";
    atomstyle = "ellipsoid";
     */
    //for(int i=0;i<3;i++) boundary.push_back("p");
    Nt=0;
    Np = 0;
}



// ---------------------------------------------------------------
// Class destructor
Particles::~Particles()
{
    //delete lmp;
}




// ---------------------------------------------------------------
// Add particle type as defined from input
void Particles::addtype(std::string newname)
{
    typenames.push_back(newname);
    Nt++;
}






// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Particles::printall()
{
	fprintf(screen,"\n---------ALL ABOUT PARTTYPE----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}