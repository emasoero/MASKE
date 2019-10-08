#include "interactions.h"
#include "particles.h"
#include "error.h"
//#include "universe.h"
//#include "DTnucleate.h"
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
*/
// for cluster Rocket
#include <string.h>


using namespace MASKE_NS;


Interactions::Interactions(MASKE *maske) : Pointers(maske)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    /*MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating lammpsIO class\n");

    lammps_active=false;
    units = "real";
    atomstyle = "ellipsoid";
     */
    Nc = 0;
}




// ---------------------------------------------------------------
// Class destructor
Interactions::~Interactions()
{
    //delete lmp;
}




// ----------------------------------------------------------------
// Builds incrementally the string to pass to lammps the pair_style
void Interactions::writestyle(std::string str)
{
    if (strlen(stylestr.c_str())==0) stylestr = str;
    else stylestr = stylestr + " " + str;
}




// ----------------------------------------------------------------
// Adds pair coefficients to the matrix of pair coefficients
void Interactions::addpcoeff(std::string p1name,std::string p2name, std::string coeffstr)
{
    // assign size of coefficient matrix (actually a vector) if this is the first coefficient
    if (Nc==0) Mcoeff.resize( (particles->Nt) * (particles->Nt) );
    
    //find id of particle and put vector coefficient string in the matrix
    int ptid1 = -1;
    int ptid2 = -1;
    
    for (int i=0; i<particles->Nt; i++) {
        if (strcmp(p1name.c_str(),(particles->typenames[i]).c_str()) == 0) ptid1 = i;
        if (strcmp(p2name.c_str(),(particles->typenames[i]).c_str()) == 0) ptid2 = i;
    }
    
    if (ptid1==-1 || ptid2==-1) {
        std::string msg = "ERROR: one of the particle types indicated for the interaction coefficients was not defined beforehand";
        error->errsimple(msg);
    }
    Mcoeff[ptid1*particles->Nt+ptid2]=coeffstr;
    Mcoeff[ptid2*particles->Nt+ptid1]=coeffstr;
    
    Nc++;
}









// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Interactions::printall()
{
	fprintf(screen,"\n---------ALL ABOUT Interactions----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
