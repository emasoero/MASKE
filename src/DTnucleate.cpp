#include "DTnucleate.h"
#include "universe.h"
#include "lammpsIO.h"
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

*/
#include <string.h>

using namespace MASKE_NS;

DTnucleate::DTnucleate(MASKE *maske) : Pointers(maske)
{
    // The inputcprs class is run by all processors in the subcommunicator. This reads the id of the processor
    //MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    //if (me == MASTER) fprintf(screen,"Generating universe class\n");

    Ntrans = 0;
    Ntsc=0;

}



// ---------------------------------------------------------------
// Class destructor
DTnucleate::~DTnucleate()
{
    
}





// ---------------------------------------------------------------
// Fucntion called at the inputcprs-reading stage, to inform all the processors about all the possible nucleation events
void DTnucleate::addtrans(std::string tname, std::string scname, std::string freqtypestr, double freqnum, double leval, double seval)
{
    Ntrans++;
    Tnames.push_back(tname);
    scnames.push_back(scname);
    if (strcmp(freqtypestr.c_str(), "step") == 0) steptimeptr.push_back(&(msk->doublestep));
    if (strcmp(freqtypestr.c_str(), "time") == 0) steptimeptr.push_back(&(msk->tempo));
    freqeval.push_back(freqnum);
    lasteval.push_back(leval);
    starteval.push_back(seval);
    
    
    
    //just a check
    /*int nprocs,subme,me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    fprintf(screen,"P%d: Nucleation \"%s\" is evaluated every %f %s\n",me,Tnames[Ntrans-1].c_str(),*steptimeptr[Ntrans-1],freqtypestr.c_str());
    */
    
}



// ---------------------------------------------------------------
// Copy subcomm-specific tranistion info into vectors, and create lammps instances as needed
void DTnucleate::addtransSC(int tn)
{
    // record to local subcommunicator
    Ntsc++;
    TnamesSC.push_back(Tnames[tn]);
    steptimeptrSC.push_back(&steptimeptr[tn]);
    freqevalSC.push_back(freqeval[tn]);
    lastevalSC.push_back(lasteval[tn]);
    startevalSC.push_back(starteval[tn]);
    
    // start a lammps instance only if needed by the current type of process/transition
    lammpsIO->create();
    
    // find the next evaluation step/time of the transition
    nextevalSC.push_back(startevalSC[Ntsc-1]);
    while ( nextevalSC[Ntsc-1] < **steptimeptrSC[Ntsc-1] ) {
        nextevalSC[Ntsc-1] += freqevalSC[Ntsc-1];
    }
    
    
    //just a check
    int nprocs,subme,me;
     MPI_Comm_rank(MPI_COMM_WORLD,&me);
     MPI_Comm_size(universe->subcomm,&nprocs);
     MPI_Comm_rank(universe->subcomm, &subme);
     fprintf(screen,"P%d: Nucleation \"%s\" allocated to rank %d in subcomm with %d procs. Its attributes are: step/time = %f , nexteval = %f \n",me,TnamesSC[Ntsc-1].c_str(),subme,nprocs,**steptimeptrSC[Ntsc-1],nextevalSC[Ntsc-1]);
    
}






// ---------------------------------------------------------------
// Checks whether the transition is currently active or not (depending on simulation step or KMC time)
bool DTnucleate::isactive(int nid)
{
    
     if (nextevalSC[nid] <= **steptimeptrSC[nid]) {
        nextevalSC[nid] += freqevalSC[nid];
        return true;
    }
    else {
        return false;
    }
}




// ---------------------------------------------------------------
// Runs the trial nucleations and updates the rate vector in class cipresso accordingly
void DTnucleate::run(int nid)
{
    //just a check
    int subme,me;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_rank(universe->subcomm, &subme);
    if (subme==0) fprintf(screen,"P%d: Running %s at step %d and time %f\n",me,TnamesSC[nid].c_str(),msk->step,msk->tempo);
    
}










    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void DTnucleate::printall()
{
	fprintf(screen,"\n---------ALL ABOUT DTnucleate----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
