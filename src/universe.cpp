#include "universe.h"
#include "error.h"
#include "randm.h"
#include "output.h"
//#include "lammpsIO.h"
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



Universe::Universe(MASKE *maske) : Pointers(maske)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    
    if (me == MASTER) fprintf(screen,"Generating universe class\n");

    SCnames.clear();
    SCnp.clear();
    nsc = 0;
    color = 0;
    key = 0;
    flampID = -1;
    flampSC = -1;

}



// ---------------------------------------------------------------
// Class destructor
Universe::~Universe()
{
    delete [] color_each;
}






// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Universe::addsubcomm(std::string name, int np, int lmpyn, int seme)
{
    //nsc ++;
    SCnames.push_back(name);
    SCnp.push_back(np);
    SClmp.push_back(lmpyn);
    SCseme.push_back(seme);
    
}





// ---------------------------------------------------------------
// Creates the universe of sub-communicators
void Universe::create()
{
    if (me == MASTER) fprintf(screen,"\nCreating the universe of sub-communicators as defined by the user in the inputcprs file (or in the restart file)...\n");
    
    // Create cumulative number of processors vector
    std::vector<int> cumnp;
    cumnp.push_back(SCnp[0]);
    for (int i=1; i<nsc; i++) cumnp.push_back(cumnp[i-1]+SCnp[i]);
    
    // Find the subcomm to which the current processor pertains. Assign it a color (subcomm number) and key (rank in the subcomm)
    for (int i=0; i<nsc; i++){
        if (SClmp[i]==1 && flampID<0 && flampSC<0){
            flampID = cumnp[i] - SCnp[i];
            flampSC = color;
        }
        if (me >= cumnp[i]) color++;
    }
    if (me<cumnp[0]) key = me;
    else key = me - cumnp[color-1];
    
    //just a check
    //fprintf(screen,"CW_P%d is now in subcomm %d with rank %d\n",me,color,key);
    
    //splitting the communicator
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &subcomm);
    
    
    // Create vector with IDs of all submasters
    subMS.push_back(0);
    for (int i=1; i<nsc; i++) {
        subMS.push_back(subMS[i-1]+SCnp[i-1]);
    }
    
    
    
    //create log file for each processor, with subcom and rank in subcom specified
    if (msk->wplog) {
        std::string fname;
        std::ostringstream ss;
        ss << me;
        fname = "p"+ss.str()+"_S";
        ss.str("");
        ss.clear();
        ss << color;
        fname = fname + ss.str()+"_k";
        ss.str("");
        ss.clear();
        ss << key;
        fname = fname + ss.str()+".plog";
        msk->plogfname=fname;
        output->createplog(fname);
    }

    if (msk->nulog_flag) {
      std::stringstream ss;
      ss << "p" << me << ".nulog";
      msk->nulog.open(ss.str());
    }

    if (msk->speclog_flag) {
      std::stringstream ss;
      ss << "p" << me << ".speclog";
      msk->speclog.open(ss.str());
    }
    
    // All processors tell every other processor their color (subcomm id) for later use (when sending subcomm-specific stuff, e.g. random number seed)
    color_each = new int[nprocs];
    color_each[me] = color;
    MPI_Allgather(MPI_IN_PLACE,1,MPI_INT,color_each,1,MPI_INT,MPI_COMM_WORLD);
    
    // Seeding random generator for each processor
    randm->seedit(SCseme[color]);
    
    
    /*int sleeptime;
    sleeptime = 2*me;
    sleep(sleeptime);
    fprintf(screen,"\n\n Processor %d, part of subcomm %d, says colors are: \n",me,color);
    for (int i=0; i<nprocs; i++) {
         fprintf(screen,"%d, ",color_each[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
    */
     
     
     
    /*
    subcID = new int[nsc];
    subcID[0] = 0;   //the MASTER
    // if a submaster, send the MASTER your processor id
    if (me > MASTER && key == 0) {
        //send "me" to master in subcID vector at position = color
        int dest = MASTER;
        MPI_Send(&subcID[color], 1, MPI_INT, dest, MPI_ANY_TAG, MPI_COMM_WORLD);
    }
    if (me=master) {
        for all subcomms > 0, so excl master
        receive each "me" at position color in subcID vector
    }
    MPI_Barrier(MPI_COMM_WORLD);
     */
    
    
    //just a check
    /*
     int nprocs,subme;
    MPI_Comm_size(subcomm,&nprocs);
    MPI_Comm_rank(subcomm, &subme);
    fprintf(screen,"CW_P%d is now part of the subcomm with %d procs. It has subrank %d\n",me,nprocs,subme);
     */
}







    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Universe::printall()
{
	fprintf(screen,"\n---------ALL ABOUT UNIVERSE----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
    //MPI_Comm_free(&subcomm);
	
}
