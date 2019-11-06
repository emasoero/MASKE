#include "error.h"

using namespace MASKE_NS;

// ---------------------------------------------------------------
// Initialize class
Error::Error(MASKE *maske) : Pointers(maske)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating error class\n");   
}




// ---------------------------------------------------------------
// Class destructor
Error::~Error()
{
    
}




// ---------------------------------------------------------------
// This produces a simple error taking care of turning off all the parallel processes
void Error::errsimple(std::string msg)
{
   // fprintf(screen,"\nOK here\n");
   //     fflush(screen);
    
    if (me==MASTER) fprintf(screen,"\n---------------------------------------------------------------\n");
    if (me==MASTER) fprintf(screen,"%s",msg.c_str());
    if (me==MASTER) fprintf(screen,"\n---------------------------------------------------------------\n");
    MPI_Finalize();
    exit(1);
}
