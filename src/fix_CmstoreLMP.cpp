#include "fix_CmstoreLMP.h"
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

Fix_CmstoreLMP::Fix_CmstoreLMP(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    //Ntrans = 0;
    //Ntsc=0;

}



// ---------------------------------------------------------------
// Class destructor
Fix_CmstoreLMP::~Fix_CmstoreLMP()
{
    
}



// ---------------------------------------------------------------
// Compute time increment of the current process.
// If static DT just read it from vector in fix.cpp. If dynamic DT, it may need to be computed somhow, hence I defined a placeholder function.
double Fix_CmstoreLMP::getDT(int pos)
{
    double Dt=0.;
    Dt = fix->Cdt[pos];
    return Dt;
}



// ---------------------------------------------------------------
//Execute store LMP fix
// Only processors in executing subcomm know the position of this local fix in the global vectors of all fix parameters, "aC...", in fix.cpp
// So executing submaster bradcasts global ID of selected fix, and then all processors can access this fix parameters from the "aC..." vectors in fix.cpp
void Fix_CmstoreLMP::execute(int pos, int subcomm)
{
    int gpos;   // global position of local fix CmstoreLMP in fix.cpp, vectors "aC..."
    //int root;   // rank of submaster of subcomm managing this Cont fix (the only one knowing gpos..)
    
    if (universe->color == subcomm && universe->key == 0) {
        //root = me;
        gpos = fix->CglobID[pos];
    }
    
    MPI_Bcast(&gpos, 1, MPI_INT, universe->subMS[subcomm], MPI_COMM_WORLD);
    
    fprintf(screen,"PROC %d: from broadcast, global ID of fix mstoreLMP to execute is %d\n",me,gpos);
    
    int sid = fix->aCsid[gpos];   // id of multi store in store.cpp
    
    if (lammpsIO->lammps_active){
        for (int i=0; i<store->MulCmd[sid].size(); i++)
            lammpsIO->lammpsdo(store->MulCmd[sid][i]);
    }
}


    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Fix_CmstoreLMP::printall()
{
	fprintf(screen,"\n---------ALL ABOUT FIX_CmstoreLMP----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
