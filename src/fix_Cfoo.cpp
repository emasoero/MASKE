#include "fix_Cfoo.h"
#include "fix.h"
//#include <sstream>
/*#include "universe.h"
#include "chemistry.h"
#include "error.h"
#include "lammpsIO.h"
#include "error.h"
#include "output.h"
#include "krun.h"
#include "solution.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

using namespace MASKE_NS;

Fix_Cfoo::Fix_Cfoo(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    //MPI_Comm_rank(MPI_COMM_WORLD, &me);
    //Ntrans = 0;
    //Ntsc=0;

}



// ---------------------------------------------------------------
// Class destructor
Fix_Cfoo::~Fix_Cfoo()
{
    
}



// ---------------------------------------------------------------
// Compute time increment of the current process
double Fix_Cfoo::getDT(int pos)
{
    double Dt=0.;
    Dt = fix->Cdt[pos];
    return Dt;
}


    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Fix_Cfoo::printall()
{
	fprintf(screen,"\n---------ALL ABOUT FIX_CFOO----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
