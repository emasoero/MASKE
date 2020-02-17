#include "fix_nufeb.h"
#include "fix.h"
#include "store.h"
#include "lammpsIO.h"
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

Fix_nufeb::Fix_nufeb(MASKE *maske) : Pointers(maske)
{
  // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
  //MPI_Comm_rank(MPI_COMM_WORLD, &me);
  //Ntrans = 0;
  //Ntsc=0;
}

// ---------------------------------------------------------------
// Class destructor
Fix_nufeb::~Fix_nufeb()
{
    
}

// ---------------------------------------------------------------
// Initialise the fix, see decription in fix_nufeb.h
void Fix_nufeb::init(int pos)
{
  int sid = fix->Csid[pos];
  for (int i=0; i<store->MulCmd[sid].size(); i++) {
    lammpsIO->lammpsdo(store->MulCmd[sid][i]);
  }
}

// ---------------------------------------------------------------
// Compute time increment of the current process
double Fix_nufeb::getDT(int pos)
{
  return fix->Cdt[pos]*fix->Csteps[pos];
}
    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Fix_nufeb::printall()
{
  fprintf(screen,"\n---------ALL ABOUT FIX_CFOO----------\n");
  //fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
  fprintf(screen,"---------------------------------------\n\n");
}

// ---------------------------------------------------------------
// Execute nufeb
void Fix_nufeb::execute(int pos)
{
  std::ostringstream ss1;
  ss1 << "timestep " << fix->Cdt[pos];
  lammpsIO->lammpsdo(ss1.str());
  std::ostringstream ss2;
  ss2 << "run " << fix->Csteps[pos];
  lammpsIO->lammpsdo(ss2.str());
}
