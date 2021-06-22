#include "setconc.h"
#include "chemistry.h"
//#include "universe.h"
//#include "lammpsIO.h"
#include "error.h"
#include "solution.h"

#include <string.h>

using namespace MASKE_NS;

// ---------------------------------------------------------------
// Initialize class
Setconc::Setconc(MASKE *maske) : Pointers(maske)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == MASTER) fprintf(screen,"Generating Setconc class\n");
}




// ---------------------------------------------------------------
// Class destructor
Setconc::~Setconc()
{
    
}





// ---------------------------------------------------------------
// record new relaxer
void Setconc::add_conc(std::string sname, int every, std::string mname,double mconc,bool flag_ctr,std::string cmol,std::string boxdV)
{
    setnames.push_back(sname);
    vevery.push_back(every);
    molnames.push_back(mname);
    molconcs.push_back(mconc);
    ctr_flags.push_back(flag_ctr);
    ctr_mols.push_back(cmol);
    vec_boxdV.push_back(boxdV);
    
    // find ID corresponding to names of molecule to be set and counterion if needed
    molID.push_back(-1);
    ctrID.push_back(-1);
    for (int i=0; i<chem->molnames.size(); i++){
        if( strcmp(mname.c_str(),chem->molnames[i].c_str())==0) molID[molID.size()-1] = i;
        if( flag_ctr && strcmp(cmol.c_str(),chem->molnames[i].c_str())==0) ctrID[molID.size()-1] = i;
    }
    if (molID.back()==-1){
        std::string msg = "ERROR: molecule "+mname+" requested in setconc command not found\n";
        error->errsimple(msg);
    }
    if (flag_ctr && ctrID.back()==-1){
        std::string msg = "ERROR: counterion "+cmol+" requested in setconc command not found\n";
        error->errsimple(msg);
    }
    
    // check that defined ion and counterion have opposite charges (or that molecule has zero charge)
    if (flag_ctr){
        // error if molecule and counterion charge have same sign
        if (chem->mol_z[molID.back()] * chem->mol_z[ctrID.back()] > 0.){
            std::string msg = "ERROR: in setconc, counterion "+cmol+" has same charge as molecule "+mname+" \n";
            error->errsimple(msg);
        }
    
        // error if molecule is charged and counterion is not
        if (chem->mol_z[molID.back()] != 0  && chem->mol_z[ctrID.back()]==0.){
            std::string msg = "ERROR: in setconc, counterion "+cmol+" is neutral, but molecule "+mname+" is charged\n";
            error->errsimple(msg);
        }
    }
}


// ---------------------------------------------------------------
// set concentration of given molecule and possibly counterion
 void Setconc::exec(int i)
{
     double c_old = chem->mol_cins[molID[i]];
     double cdV_old = chem->mol_cindV[molID[i]];
     double n_old = chem->mol_nins[molID[i]];
     double ndV_old = chem->mol_nindV[molID[i]];
     
     chem->mol_cins[molID[i]] = molconcs[i];
     chem->mol_nins[molID[i]] = chem->mol_cins[molID[i]] * solution->SVol;
     
     if (strcmp(vec_boxdV[i].c_str(),"box+dV")==0) {
         chem->mol_cindV[molID[i]] = molconcs[i];
         chem->mol_nindV[molID[i]] = chem->mol_cindV[molID[i]] * solution->dVSVol;
     }
     
     // same for counterion
     if (ctr_flags[i]){
         double dn;
         dn =  (chem->mol_nins[molID[i]] - n_old) * (-chem->mol_z[molID[i]]/chem->mol_z[ctrID[i]]);
         chem->mol_cins[ctrID[i]] += dn/solution->SVol;
         chem->mol_nins[ctrID[i]] += dn;
         
         if (strcmp(vec_boxdV[i].c_str(),"box+dV")==0) {
             dn =  (chem->mol_nindV[molID[i]] - ndV_old) * (-chem->mol_z[molID[i]]/chem->mol_z[ctrID[i]]);
             chem->mol_cindV[ctrID[i]] += dn/solution->dVSVol;
             chem->mol_nindV[ctrID[i]] += dn;
         }
     }
    
}


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Setconc::printall()
{
    fprintf(screen,"\n---------ALL ABOUT SETCONC----------\n");
    //fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
    fprintf(screen,"---------------------------------------\n\n");
    
}
