#include "lammpsIO.h"
#include "universe.h"
#include "simbox.h"
#include "interactions.h"
#include "particles.h"
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


LammpsIO::LammpsIO(MASKE *maske) : Pointers(maske)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating lammpsIO class\n");

    lammps_active=false;
    units = "real";
    atomstyle = "ellipsoid";
    timestep = 1.;
    lmpThSt = "thermo_style custom step atoms";
}



// ---------------------------------------------------------------
// Class destructor
LammpsIO::~LammpsIO()
{
    delete lmp;
}



// ----------------------------------------------------------------
// Create a new Lammps instance
void LammpsIO::create()
{
    // if lammps has not been already activated by this proc (just for safety) and if the processor pertains to a subcomm for which lammps was turend on in the input file..
    if (!lammps_active && universe->SClmp[universe->color]==1) {
        
        // option to print or not print lammps output to screen (useful for debugging but littering MASKE's output otherwsie)
        bool noscreen = true;
        int argn = 0;
        char **words;
        if (noscreen) {
            words = (char**)calloc (3,sizeof(char*));
            words[0] = "placeholder";
            words[1] = "-screen";
            words[2] = "none";
            argn = 3;
        }
        else words = NULL;

        lmp = new LAMMPS_NS::LAMMPS(argn,words,universe->subcomm);
        //lmp = new LAMMPS_NS::LAMMPS(0,NULL,universe->subcomm);
        lammps_active=true;
        
        delete [] words;
       
        
        fprintf(screen,"Proc[%d]: Lammps instance created, as part of subcomm %s\n",me,universe->SCnames[universe->color].c_str());
        
        // creates subcomm-specific lammps logfile (and temporary dumpfile name, currently turned off)
        std::string todo;
        todo = "log log."+universe->SCnames[universe->color]+".lammps";
        lammpsIO->lammpsdo(todo);
        //tdump_fname = universe->SCnames[universe->color]+".tdump";
        
        
    /*
        std::string instr = "units "+units;
        fprintf(screen,"\n Proc[%d]: units command is: %s\n",me,instr.c_str());
        lmp->input->one(instr.c_str());
        
        instr = "atom_style "+atomstyle;
        //fprintf(screen,"\n atom_style command is: %s\n",instr.c_str());
        lmp->input->one(instr.c_str());
        
        instr = "boundary "+simbox->boundary[0]+" "+simbox->boundary[1]+" "+simbox->boundary[2];
        //fprintf(screen,"\n Boundary command is: %s\n",instr.c_str());
        lmp->input->one(instr.c_str());
        
        instr = "pair_style "+interact->stylestr;
        //fprintf(screen,"\n Pair_style command is: %s\n",instr.c_str());
        lmp->input->one(instr.c_str());
        
        ss<< "timestep "<<timestep;
        instr = ss.str();
        ss.str("");
        //fprintf(screen,"\n Proc[%d]: timestep command is: %s\n",me,instr.c_str());
        lmp->input->one(instr.c_str());
        
        //lammps creates the simulation box
        ss<<"region simboxENRICOswwwlp prism "<<simbox->xyzlo[0]<<" "<<simbox->xyzhi[0]<<" "<<simbox->xyzlo[1]<<" "<<simbox->xyzhi[1]<<" "<<simbox->xyzlo[2]<<" "<<simbox->xyzhi[2]<<" "<<simbox->xyztri[0]<<" "<<simbox->xyztri[1]<<" "<<simbox->xyztri[2];
        instr = ss.str();
        ss.str("");
        fprintf(screen,"\n Proc[%d]: region command is: %s\n",me,instr.c_str());
        lmp->input->one(instr.c_str());
        ss<<"create_box "<<particles->Nt<<" simboxENRICOswwwlp";
        instr = ss.str();
        ss.str("");
        fprintf(screen,"\n Proc[%d]: create_box command is: %s\n",me,instr.c_str());
        lmp->input->one(instr.c_str());
        
        
        // pair coefficient passed to lammps
        for (int i=0; i<particles->Nt; i++) {
            for (int j=i; j<particles->Nt; j++) {
                ss<<"pair_coeff "<<i+1<<" "<<j+1<<" "<<interact->Mcoeff[i*particles->Nt+j];
                instr = ss.str();
                ss.str("");
                if (me == MASTER) fprintf(screen,"\n pair_coeff command is: %s\n",instr.c_str());
                lmp->input->one(instr.c_str());
            }
        }
        
        */
        
    }
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void LammpsIO::lammpsdo(std::string todo)
{
   lmp->input->one(todo.c_str());
    
}






// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void LammpsIO::printall()
{
	fprintf(screen,"\n---------ALL ABOUT LAMMPSIO----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
