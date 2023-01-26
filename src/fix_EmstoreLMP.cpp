#include "fix_EmstoreLMP.h"
#include "fix.h"
#include "lammpsIO.h"
#include "universe.h"
#include "store.h"
//#include <sstream>
/*
#include "chemistry.h"
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

Fix_EmstoreLMP::Fix_EmstoreLMP(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    //Ntrans = 0;
    //Ntsc=0;

}



// ---------------------------------------------------------------
// Class destructor
Fix_EmstoreLMP::~Fix_EmstoreLMP()
{
    
}




// ---------------------------------------------------------------
//Execute store LMP fix
// Only processors in executing subcomm know the position of this local fix in the global vectors of all fix parameters, "aC...", in fix.cpp
// So executing submaster bradcasts global ID of selected fix, and then all processors can access this fix parameters from the "aC..." vectors in fix.cpp
void Fix_EmstoreLMP::execute(int gpos)
{
    // gpos = global position of local fix EmstoreLMP in fix.cpp, vectors "aE..."
    
    int sid = fix->aEsid[gpos];   // id of multi store in store.cpp
    
    if (lammpsIO->lammps_active){
        for (int i=0; i<store->MulCmd[sid].size(); i++)
            lammpsIO->lammpsdo(store->MulCmd[sid][i]);
    }
}


    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Fix_EmstoreLMP::printall()
{
	fprintf(screen,"\n---------ALL ABOUT FIX_EmstoreLMP----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
