#include "maske.h"

/*
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
*/

#include "memory.h"
#include "error.h"
#include "inputmsk.h"
#include "universe.h"
#include "DTnucleate.h"
#include "lammpsIO.h"
#include "simbox.h"
#include "particles.h"
#include "chemistry.h"
#include "interactions.h"
#include "solution.h"
#include "fix.h"
#include "fix_delete.h"
#include "fix_Cfoo.h"
#include "relax.h"
#include "krun.h"
#include "randm.h"
#include "output.h"
#include "fix_nucleate.h"
#include "block.h"
#include "store.h"

#ifdef MASKE_WITH_NUFEB
#include "fix_nufeb.h"
#endif

using namespace MASKE_NS;

// ---------------------------------------------------------------
// Initialize class
MASKE::MASKE(int narg, char **arg)
{
    screen = NULL;
    plog = NULL;
    screen=stdout;
    wplog = false;
    nulog_flag = false;

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me==MASTER) {
        fprintf(screen,"\n---------STARTING SIMULATION----------\n");
        fprintf(screen,"argc = %d\n",narg);
        for (int i=0; i<narg; i++) {
            fprintf(screen,"arg[%d] = %s\n",i,arg[i]);
        }
        fprintf(screen,"---------------------------------------\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    Nsteps = 0;
    step = 0;
    tempo = 0;
    doublestep = (double)step;
    kB = 1.;
    hpl = 1.;
    
    Rtypes.clear();
    Ttypes.clear();
	
	memory = new Memory(this);
    error = new Error(this);
	inputmsk = new Inputmsk(this,narg,arg);
    universe = new Universe(this);
    nucleate = new DTnucleate(this);
    lammpsIO = new LammpsIO(this);
    simbox = new Simbox(this);
    particles = new Particles(this);
    chem = new Chemistry(this);
    interact = new Interactions(this);
    solution = new Solution(this);
    fix = new Fix(this);
    fix_del = new Fix_delete(this);
#ifdef MASKE_WITH_NUFEB
    fix_nufeb = new Fix_nufeb(this);
#endif
    krun = new Krun(this);
    randm = new Randm(this);
    output = new Output(this);
    fix_cfoo = new Fix_Cfoo(this);
    relax = new Relax(this);
    fix_nucl = new Fix_nucleate(this);
    block = new Block(this);
    store = new Store(this);
}

// ---------------------------------------------------------------
// Class destructor
MASKE::~MASKE()
{
    if (me==MASTER) fprintf(screen,"Deleting maske class\n");
    MPI_Barrier(MPI_COMM_WORLD);
    if (me==MASTER) {
        fprintf(screen,"Deleting inputmsk class\n");
        inputmsk->printall();
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete inputmsk;
    if (me==MASTER) fprintf(screen,"Deleting error class\n");
    MPI_Barrier(MPI_COMM_WORLD);

    delete error;
    if (me==MASTER) {
        fprintf(screen,"Deleting DTnucleate class\n");
        nucleate->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete nucleate;
    if (me==MASTER) {
        fprintf(screen,"Deleting lammpsIO class\n");
        lammpsIO->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete lammpsIO;
    if (me==MASTER) {
        fprintf(screen,"Deleting chemistry class\n");
        chem->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete chem;
    if (me==MASTER) {
        fprintf(screen,"Deleting particles class\n");
        particles->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete particles;
    if (me==MASTER) {
        fprintf(screen,"Deleting simbox class\n");
        simbox->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete simbox;
    if (me==MASTER) {
        fprintf(screen,"Deleting interactions class\n");
        interact->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete interact;
    if (me==MASTER) {
        fprintf(screen,"Deleting solution class\n");
        solution->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete solution;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix class\n");
        fix->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix;
    if (me==MASTER) {
        fprintf(screen,"Deleting krun class\n");
        krun->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete krun;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix_delete class\n");
        fix_del->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix_del;
    if (me==MASTER) {
        fprintf(screen,"Deleting randm class\n");
        randm->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete randm;
    if (me==MASTER) {
        fprintf(screen,"Deleting Cfoo class\n");
        fix_cfoo->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix_cfoo;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix_cfoo class\n");
        relax->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef MASKE_WITH_NUFEB
    delete fix_nufeb;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix_nufeb class\n");
        relax->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    delete relax;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix_nucleate class\n");
        fix_nucl->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix_nucl;
    if (me==MASTER) {
        fprintf(screen,"Deleting block class\n");
        block->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete block;
    if (me==MASTER) {
        fprintf(screen,"Deleting store class\n");
        store->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete store;
    
    if (me==MASTER) {
        fprintf(screen,"Deleting universe class\n");
        universe->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    /*fflush(screen);
    int sleeptime = 2;
    sleep(sleeptime);
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(screen,"\n IN0 PROC %d: eddig jo",me);
    sleep(sleeptime);
    MPI_Barrier(MPI_COMM_WORLD);
    fflush(screen);*/
    MPI_Comm_free(&(universe->subcomm));

    //sleep(sleeptime);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    
      delete universe;
    
    //MPI_Comm_free(&(universe->subcomm));

    //sleeptime = 2;
    //sleep(sleeptime);
   /* MPI_Barrier(MPI_COMM_WORLD);
    fprintf(screen,"\n IN PROC %d: eddig jo",me);
    fflush(screen);
    fflush(screen);
    sleep(4);
    fflush(screen);
    MPI_Barrier(MPI_COMM_WORLD);
    fflush(screen);
*/

	
}


// ---------------------------------------------------------------
// The main loop    //OBSOLETE
void MASKE::mainloop()
{
    
    
    // Safety check, run by all processors in the universe
    if (Nsteps <= 0) {
        std::string msg = "ERROR: the number of fundamental steps \"Nsteps\" for the main loop is not a positive integer. \n Did you remember to set a positive value of this number in the inputcprs file?";
        error->errsimple(msg);
    }
    
    //just a check
    /*
     if (me==MASTER) {
        for (int i=0; i<particles->Nt; i++) {
            for (int j=0; j<particles->Nt; j++) {
                fprintf(screen,"Coeff string for %s and %s is %s \n",(particles->typenames[i]).c_str(),(particles->typenames[j]).c_str(),(interact->Mcoeff[i*particles->Nt+j]).c_str());
            }
        }
    }
    */

    
    
    
    // If this is not a restarted simulation, save the very firt state of the system
    if (inputmsk->isrestart==false) {
        //save system at step0 and corresponging tempo
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Run the main loop (the initial value of "step" is not defined here because it might be read from a restart
    while (step<Nsteps) {
        
        doublestep = (double) step;
        
        if (me==MASTER) fprintf(screen,"\nHere we are: step %d\n",step);
        
        // Each sub-communicator runs all its discrete kinetic transitions. Each of these actions pushes back rates and quantities associated with the transition (e.g. particle position and shape in case of nucleation, or new state in case of creep transition). Each command below will correpsond to an object where the needed info where preemtly stored
        
        // Run active nucleation events now...
        for (int nid=0; nid < nucleate->Ntsc; nid++) {
            if (nucleate->isactive(nid)) nucleate->run(nid);
            // the code requires a lammps-dump-like xyz file with additional columns for the molecular composition of the particles. The xyz is used to generate the particles in cipresso and to create a properly formatted data file that is then loaded in lammps (this is needed to preserve the ids, which are otherwise lost if the particles are created in lammps instead of loaded from the data file)
            // reset timestep in lammps to match the one read from the xyz file
            // printing xyz file at regular intervals
            // implement the possibility to define n lattices with specific names or/and ids, so that lammps can invoke nucleations that are based on different lattices if needed
            // create a background object that can contain a vector of background fields
            // import and write the background
            // within nucleate, add arguments to freeze particle types and to manage neighbours
        }
        /*
        for (int did=0; did < deletepart->Ntsc; did++) {
            if (creepOSC->isactive(ncosc)) creepOSC->run(ncosc);
        }
        for (int ncqs=0; ncqs < creepQS->Ntrans; ncqs++) {
            if (creepQS->isactive(ncqs)) creepQS->run(ncqs);
        }
        for (int ncosc=0; ncosc < DTcreepOSC->Ntrans; ncosc++) {
            if (creepOSC->isactive(ncosc)) creepOSC->run(ncosc);
        }*/
        MPI_Barrier(MPI_COMM_WORLD);
        
        
        
        // The MASTER collects all the possible transition rates and some sort of pointers that allow to go back to the chosen transition in the chosen processor
        
        
        
        
        // The MASTER chooses the winning transition and updates the system and time accordingly. If needed, saves the foreground
        tempo += 0.055;
        step++;
        
        
        // Within a loop of incremental time steps, each processor runs all the continuum processes that are active at the current step, or active between the previous and just-updated time
    }
}

// ---------------------------------------------------------------
// Print stuff
void MASKE::printall()
{
    if (me==MASTER) {
        fprintf(screen,"\n---------END OF SIMULATION----------\n");
        //fprintf(screen,"water to cement ratio =  %f\n",wc);
        fprintf(screen,"---------------------------------------\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
