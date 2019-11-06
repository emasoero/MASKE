#include "relax.h"
//#include "chemistry.h"
//#include "universe.h"
#include "lammpsIO.h"
//#include "error.h"

#include <string.h>

using namespace MASKE_NS;

// ---------------------------------------------------------------
// Initialize class
Relax::Relax(MASKE *maske) : Pointers(maske)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == MASTER) fprintf(screen,"Generating fix_Relax class\n");
}




// ---------------------------------------------------------------
// Class destructor
Relax::~Relax()
{
    
}





// ---------------------------------------------------------------
// record new relaxer
void Relax::add_rlx(std::string rid,int every,std::string relaxer,std::string rstring, std::string rstyle, std::string rmodify)
{
    rlxID.push_back(rid);
    rlx_every.push_back(every);
    rlx_string.push_back(rstring);
    rlx_relaxers.push_back(relaxer);
    if (strcmp(relaxer.c_str(),"minimize")==0) {
        // only minimization requires min_style and min_modify
        rlx_style.push_back(rstyle);
        rlx_modify.push_back(rmodify);
    }
    else {
        rlx_style.push_back("");
        rlx_modify.push_back("");
    }
 
}


// ---------------------------------------------------------------
// write new entry in dump file
 void Relax::dorelax(int i)
{
    std::string tolmp;
    
    //if relaxer is a minimze, preamble required to set style and modify
    if (strcmp(rlx_relaxers[i].c_str(),"minimize")==0) {
        tolmp = "min_style "+rlx_style[i];
        lammpsIO->lammpsdo(tolmp);
        lammpsIO->lammpsdo(rlx_modify[i]);
    }
    
    // run minimizer
    lammpsIO->lammpsdo(rlx_string[i]);
    
}


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Relax::printall()
{
    fprintf(screen,"\n---------ALL ABOUT FIX_RELAX----------\n");
    //fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
    fprintf(screen,"---------------------------------------\n\n");
    
}
