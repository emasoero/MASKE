#include "fix_delete.h"
#include "fix.h"
//#include <sstream>
#include "universe.h"
#include "chemistry.h"
#include "error.h"
#include "lammpsIO.h"
#include "error.h"
#include "output.h"
#include "krun.h"
#include "solution.h"
#include <math.h>       /* for the pow function used here
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#include <string.h>

using namespace MASKE_NS;

Fix_delete::Fix_delete(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    EVpID = -1;
    SAR = nullptr;
    //Ntrans = 0;
    //Ntsc=0;
}



// ---------------------------------------------------------------
// Class destructor
Fix_delete::~Fix_delete()
{
    
}





// ---------------------------------------------------------------
// Initialise the fix, here it simply counts the number of particles to be deleted by the current fix. In the future I might implement a region as well, so that deletion will be for particles of selected tpye in given region block (region blocks can dynamically follow the simulation box deformations)
void Fix_delete::init(int pos)
{
    //number of delete events required, then recorded into fix.cpp
    std::ostringstream typ;
    typ << fix->fKMCptype[pos];
    std::string tolmp = "group gtemp type "+typ.str();
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable tempNgroup equal count(gtemp)";
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");
    double tNg = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"tempNgroup",0));
    int Ng = (int) tNg;
    lammpsIO->lammpsdo("group gtemp delete");
    lammpsIO->lammpsdo("variable tempNgroup delete");
    // the subsequent two lines will be used in fix_del to know where to start and where to finish scanning events in the assembled vectors of ids, rates, etc..
    fix->fKMCnevents[pos] = Ng;
    
    
}




// ---------------------------------------------------------------
// Sampling all delete transitions correspondin to fix in position "pos" within this subcommunicator list, recorded in the fix object
void Fix_delete::sample(int pos)
{
    
   
    
    // USEFUL SHORTCUTS for current processor key and number of procs in local subcomm
    key = universe->key;
    nploc =universe->SCnp[universe->color];
    
    mid = fix->fKMCmid[pos];    // mechanism ID
    
    // screen output to debug
    int sleeptime = me;
    if (msk->wplog) {
        std::string msg = "\nTEMPO ";
        std::ostringstream ss;
        ss << msk->tempo;
        msg = msg+ss.str()+": sampling fix_delete \""+fix->fKMCname[pos]+"\"";
        output->toplog(msg);
    }
    
    
    
    
    
    
    // EXTRACTING IDs, RADIUS, AND ENERGIES FROM LAMMPS
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    std::string tolmp;
    tolmp = "compute tempID all property/atom id"; // temp compute to get id of all atoms
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "compute tempRAD all property/atom radius"; // temp compute to get atom radii 
    lammpsIO->lammpsdo(tolmp);
    
    std::ostringstream typ;
    typ << fix->fKMCptype[pos];
    tolmp = "group gtemp type "+typ.str();
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable tempNgroup equal count(gtemp)";
    //tolmp = "variable tempNgroup equal count("+fix->fKMCgroup[pos]+")";  // temp variable counting atoms of type in fix
    lammpsIO->lammpsdo(tolmp);
    
    // NB: IF POTENTIAL IS NOT PAIRWISE, but e.g. three body, will need to change relationship between pe/atom extracted here and DU computed later for the rate. If potential is a mix of pair and other terms, e.g. a generic one for molecular interactions, I will need to compute energies per particle or molecule in a different way, by actually deleting it, measuring DUtot, and then putting it back - more expensive
    tolmp = "compute tempPE all pe/atom";  // temp compute reading energy per atom
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "compute tempType all property/atom type";  // temp compute reading type per atom
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "dump tdID all custom 1 dump.temp_"+universe->SCnames[universe->color]+" c_tempID c_tempRAD c_tempPE type x y z";    //temp dump to update variables and computes
    lammpsIO->lammpsdo(tolmp);
    
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            tolmp = "compute tempPAT all property/local patom1 patom2 ptype1 ptype2";
            lammpsIO->lammpsdo(tolmp);
            tolmp = "compute tempDIST all pair/local dist";
            lammpsIO->lammpsdo(tolmp);
            tolmp = " dump tdLID all local 1 dump.tempL_"+universe->SCnames[universe->color]+" index c_tempPAT[1] c_tempPAT[2] c_tempPAT[3] c_tempPAT[4] c_tempDIST";
            lammpsIO->lammpsdo(tolmp);
        }
    }
    
    lammpsIO->lammpsdo("run 1");     // a run1 in lammps to dump the temp and so prepare variables and computes

    
    
    
    // extracting unsorted group atom ids and energies
    // NB: lammps computes of properties per atom per group store the results in arrays whose size is natoms (number of atoms in the simulations, including those not in the group). Entries associated to atoms not in the group are set to 0, which I spot and get rid of
    double tNg = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"tempNgroup",0));
    Ng = (int) tNg;
    lammpsIO->lammpsdo("group gtemp delete");
    lammpsIO->lammpsdo("variable tempNgroup delete");

    
    
    //safety check: if not just_reset and Ng not equal to nevents, error
    if (Ng != fix->fKMCnevents[pos]) {
        std::ostringstream wd1,wd2,wd3,wd4,wd5;
        wd1<<me; wd2<<universe->color; wd3<<msk->tempo; wd4<<Ng; wd5<<fix->fKMCnevents[pos];
        std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: fix_delete.cpp found "+wd4.str()+" atoms in fix "+fix->fKMCname[pos]+", whereas krun.cpp previously found "+wd5.str()+" atoms in same fix. This should never happen. Something is wrong with the source code... \n****************************************************\n***********************************************\n";
        fprintf(screen,"%s",msg.c_str());
        //error->errsimple(msg);
    }
    
    natoms = static_cast<int> (lammpsIO->lmp->atom->natoms); // total number of atoms in lammps, all groups all fixes)
    nlocal = static_cast<int> (lammpsIO->lmp->atom->nlocal); // number of atoms in current processore, all types all fixes)
    
    //double *aID = NULL;
    //double *aID = (double*) malloc(natoms * sizeof(double));
    //aID = new double[natoms];     // vector of unsorted IDs
    aID = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempID",1,1));
      // NB: each processor in subcom pulls out the ID of their atoms. We will put them all into a single vector, IDuns, to be managed by the submaster. The way that seems to work is to scan aID of each processor looking for the first nlocal atoms with non-zero id. Ids after those are random. The first nonzero nlocal ids are passed to the submaster, which eventually sorts them.
    //int *aaID;      // array with size natoms, with IDs from lammps for ALL natoms
    
    //aaID = new int[natoms];  // the array containing all IDs in lammps, irrespective of group
    //int *aaID = ((int *) lammps_extract_atom(lammpsIO->lmp,(char *)"id"));
    
    //aR = new double[natoms];     // arrays of radii and energies
    aR = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempRAD",1,1));
    
    //aE = new double[natoms];
    aE = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempPE",1,1));
    
    //atype = new int[natoms];
    //int *atype = ((int *) lammps_extract_atom(lammpsIO->lmp,(char *)"type"));
    //atype = new double[natoms];
    atype = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempType",1,1));
    
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            // PAT is a local compute (style = 2 in lammps library) of array (i.e. matrix) type (type = 2 in lammps library)
            locLMP = ((double **) lammps_extract_compute(lammpsIO->lmp,(char *) "tempPAT",2,2));
            
            // DIST is a local compute (style = 2 in lammps library) of vector type (type = 1 in lammps library)
            aDIST = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempDIST",2,1));
            
            nlocR = *((int*)lammps_extract_compute(lammpsIO->lmp,(char *) "tempPAT",2,0));
            
        }
    }
    
    
    /*if (me==3){
        fprintf(screen,"\n Proc %d: \n",me);
        for (int j=0; j<natoms; j++){
            //fprintf(screen,"%d %d %d %f %f\n",j,(int)aID[j],atype[j],aR[j],aE[j]);
            //fprintf(screen,"%d %d %f %f\n",j,(int)aID[j],aR[j],aE[j]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);*/
    
    
    //printing group-specific aIDs and all IDs in lammps vectors, for debugging
    if (msk->wplog) {
        std::string msg = "Total number of atoms in LAMMPS: ";
        std::ostringstream ss;
        ss << natoms;       msg = msg+ss.str(); ss.str("");   ss.clear();
        msg = msg+"\nNumber of atoms in current processor: ";
        ss << nlocal;       msg = msg+ss.str(); ss.str("");   ss.clear();
        output->toplog(msg);
        msg="";
        output->toplog("\npos aiD aR aE");
        {
            int naf = 0;  // counter of number of atoms in current processor
            int i = 0;
            while (naf < nlocal){
                if (aID[i]>0){
                    naf++;
                    if ((int)atype[i]==fix->fKMCptype[pos]){
                        ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
                        ss << aID[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << atype[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        //ss << aaID[i];  msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << aR[i];    msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << aE[i];    msg += ss.str();      ss.str("");   ss.clear();
                        output->toplog(msg);
                        msg="";
                    }
                }
                i++;
            }
        }
        
        /*output->toplog("\nALL IDs FROM LAMMPS \npos aiD atype aR aE");
        for (int i=0; i<natoms; i++){
            ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
            ss << aID[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << atype[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            //ss << aaID[i];  msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << aR[i];    msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << aE[i];    msg += ss.str();      ss.str("");   ss.clear();
            output->toplog(msg);
            msg="";
        }*/
        
        if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
            if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
                msg = msg+"\nNumber of interacting pairs in current processor: ";
                ss << nlocR;       msg = msg+ss.str(); ss.str("");   ss.clear();
                output->toplog(msg);
                msg="";
                output->toplog("\npos id1 id2 type1 type2 dist");
                {
                    for(int i=0; i<nlocR; i++){
                        ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
                        ss << locLMP[i][0];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << locLMP[i][1];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << locLMP[i][2];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << locLMP[i][3];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << aDIST[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        output->toplog(msg);
                        msg="";
                    }
                }
            }
        }
       
        
        
    }
    
    //fprintf(screen,"\n Proc %d, EDDIG JO0? \n",me);
    //MPI_Barrier(MPI_COMM_WORLD);
    //sleep(10);
    
    tolmp = "undump tdID";     // removing temporary dump
    lammpsIO->lammpsdo(tolmp);
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            tolmp = "undump tdLID";     // removing temporary "local" dump for per-pair interactions
            lammpsIO->lammpsdo(tolmp);
        }
    }
    
    // END OF READING FROM LAMMPS
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    
    
    

    // RECORDING IDs, RADIUS, AND ENERGIES IN EACH PROCESSOR (cleaning out those not pertaining to current processor)
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    tID.clear();    tR.clear();    tE.clear();
    naP = 0;  // number of atoms in current processor of type invoked by current fix
    int naf = 0;  // counter of number of atoms in current processor
    {
        int i = 0;
        while (naf < nlocal){
            if (aID[i]>0){
                naf++;
                if ((int)atype[i]==fix->fKMCptype[pos]){
                    tID.push_back((int)aID[i]);
                    tR.push_back(aR[i]);
                    tE.push_back(aE[i]);
                    naP++;
                }
            }
            i++;
        }
    }

    tIDarr = new int[naP];     // copying local ID vector to array in order to communicate it to submaster later on
    for (int i =0; i<naP; i++) tIDarr[i] = tID[i];

    
    //printing proc-specific tIDs (both vecto and array), and all IDs in lammps vectors, for debugging
    if (msk->wplog) {
        std::string msg = "\nNumber of atoms to delete in this processor: ";
        std::ostringstream ss;    ss << naP;   msg = msg+ss.str();
        output->toplog(msg);
        msg="";   ss.str("");   ss.clear();
        output->toplog("\npos tiD tiDarr tR tE");
        for (int i=0; i<naP; i++){
            ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
            ss << tID[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << tIDarr[i];  msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << tR[i];    msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << tE[i];    msg += ss.str();      ss.str("");   ss.clear();
            output->toplog(msg);
            msg="";
        }
    }
    
    
    
    // for micro mechanism, create another array with IDs adn radii of ALL real particles, to be used to compute coverage fractions
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            
            tIDar.clear();    tRar.clear();
            
            nR = 0;  //number of real particles in current processor

            for (int i=0 ; i<nlocal; i++) {  //NB: local atoms may also be trial types when nucleation fixes are used
                bool flag_real = false;
                for (int j=0; j<msk->Rtypes.size(); j++) {
                    if ((int)atype[i]==msk->Rtypes[j]) flag_real = true;
                }
                if (flag_real) {
                    tIDar.push_back((int)aID[i]);
                    tRar.push_back(aR[i]);
                    nR++;
                }
            }
            
            IDar = new int[nR];
            Rar = new double[nR];
            
            for (int i =0; i<nR; i++){
                IDar[i] = tIDar[i];
                Rar[i] = tRar[i];
            }
            
            
            // printing to log file for debugging
            if (msk->wplog) {
                std::string msg = "\nNumber of real atoms in this processor (for coverage function in fix_delete.cpp): ";
                std::ostringstream ss;    ss << nR;   msg = msg+ss.str();
                output->toplog(msg);
                msg="";   ss.str("");   ss.clear();
                output->toplog("\npos IDar Rar");
                for (int i=0; i<nR; i++){
                    ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
                    ss << IDar[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                    ss << Rar[i];    msg += ss.str()+" ";  ss.str("");   ss.clear();
                    output->toplog(msg);
                    msg="";
                }
            }
        }
    }
    
    
    // END OF RECORDING LAMMPS VALUES INTO VECTORS
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    
    
    
    
    // SUBMASTER ASSEMBLES THE UNSORTED ID ARRAYS FROM EACH PROCESSOR IN AN UNSORTED VECTOR
    ids_to_submaster(pos);
    
    // FOR MICRO-PAIR MECHANISM ONLY Communicate local arrays (from neighbour list) to submaster
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            pair_arr_to_submaster(pos);
        }
    }
    
    // SUBMASTER SORTS IDs (only if fix was just reset in krun.cpp, otherwise sorted vector already exists)
    if (key==0 && krun->fMC_justreset) submaster_sort_IDs(pos);
    
    // EACH PROCESSOR MAPS THEIR ID TO THE CORRESPONDING POSITION IN THE SUBMASTER IDsrt ARRAY
    if(key==0) submaster_map_ID(pos);
    
    // FOR MICRO-PAIR MECHANISM ONLY, SUBMASTER COMPUTES COVERAGE AREAS AND PASSES THEM BACK TO SLAVES
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            if(key==0){
                // submaster computes coverage area fraction of each particle in the fix (stored in unsorted array in submaster)
                submaster_comp_cover(pos);
            }
            // submaster communicates back its areas to each slave processor in chunks of same size as local IDs
            cover_from_submaster(pos);
        }
    }
    
    
    // EACH PROCESSOR COMPUTES RATES OF EACH EVENT AND CUMULATIVES DEPENDING ON MECHANISM
    rate_each = new double[nID_each[key]];  //array with processor-specific rates. These will be passed to the submaster which will sort them to match ids into the rates vector. Then this array can be forgotten
    if (strcmp((chem->mechstyle[mid]).c_str(),"allpar")==0) {
        comp_rates_allpar(pos);
    }
    else if (strcmp((chem->mechstyle[mid]).c_str(),"allser")==0) {
        comp_rates_allser(pos);
    }
    else if (strcmp((chem->mechstyle[mid]).c_str(),"micro")==0) {
        comp_rates_micro(pos);
    }
    //else if (strcmp((chem->mechstyle[mid]).c_str(),"shvab")==0) {
      //  comp_rates_shvab(pos);
    //}
    else {
        fprintf(screen,"\n\n\n\n  ******************* \n ************************** ERROR : Mechanism style not recognized. \n ******************************* \n ******************** \n\n\n\n\n");
    }
    
    // PASS RATES TO SUBMASTER
    rates_to_submaster(pos);
    
    

    tolmp = "uncompute tempID";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "uncompute tempRAD";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "uncompute tempPE";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "uncompute tempType";
    lammpsIO->lammpsdo(tolmp);
    
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            tolmp = "uncompute tempPAT";
            lammpsIO->lammpsdo(tolmp);
            tolmp = "uncompute tempDIST";
            lammpsIO->lammpsdo(tolmp);
        }
    }
    
    delete [] rate_each;
    delete [] tIDarr;
    delete [] IDuns;
    delete [] IDpos;
    delete [] nID_each;


    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            
            if (!(SAR == nullptr)){
                free(SAR[0]);
                free(SAR);
                SAR = nullptr;
            }
            delete [] tCF;
            delete [] tGM;
            
            delete [] nlocR_each;
            delete [] SARpos;
            delete [] Dsub;
            delete [] CFuns;
            delete [] GMuns;
            delete [] fGMuns;
            
            delete [] IDar;
            delete [] Rar;
            delete [] nIDar_each;
            delete [] IDarpos;
            delete [] IDaruns;
            delete [] Raruns;

        }
    }

}




// ---------------------------------------------------------------
// assembles all IDs from each processor into block-wise array in submaster, where each block is the tIDarr in a processor
void Fix_delete::ids_to_submaster(int pos)
{
    //MPI send size of each tID to submaster. Submaster creates array with enough space and assigns positions to accept IDs.  All procs then send tIDarr to sub master
    nID_each = new int[nploc];
    nID_each[key] = naP;   //each processor in subcomm records its number of atoms (naP) at location = key of its nID_each
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            nIDar_each = new int[nploc];
            nIDar_each[key] = nR;
        }
    }
    
    if (key>0) {
        int dest = 0;
        MPI_Send(&nID_each[key], 1, MPI_INT, dest, 1, (universe->subcomm));
    }
    if (key==0 && nploc>1) {
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&nID_each[source], 1, MPI_INT, source, 1, (universe->subcomm), &status);
        }
    }
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            if (key>0) {
                int dest = 0;
                MPI_Send(&nIDar_each[key], 1, MPI_INT, dest, 1, (universe->subcomm));
            }
            if (key==0 && nploc>1) {
                for (int source=1; source<nploc; source++) {
                    MPI_Recv(&nIDar_each[source], 1, MPI_INT, source, 1, (universe->subcomm), &status);
                }
            }
        }
    }
    
    
    // check that tot number of IDs from all processors equals number of atoms in group
    if (key == 0) {
        int totIDs = 0;
        for (int i=0; i<nploc; i++) {
            totIDs = totIDs + nID_each[i];
        }
        if (totIDs != Ng) {
            std::ostringstream wd1,wd2,wd3,wd4,wd5;
            wd1<<me; wd2<<universe->color; wd3<<msk->tempo; wd4<<totIDs; wd5<<Ng;
            std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: number of IDs passed to local master "+wd4.str()+" inconsistent with number of atoms in group "+ wd5.str() +". This should never happen. Problem with the source code. \n****************************************************\n***********************************************\n";
            fprintf(screen,"%s",msg.c_str());
            // error->errsimple(msg);
        }
    }
    
    // pass local ID arrays to submaster
    IDpos = new int[nploc];  // position of local tID array in submaster's unsorted list of IDs
    IDarpos = new int[nploc];
    
    if (key==0) {
        IDpos[0]=0;
        for (int i=1; i<nploc; i++) {
            //fprintf(screen,"i =  %d  --> nID_each = %d \n",i,nID_each[i]);
            IDpos[i] =IDpos[i-1]+nID_each[i-1];
        }
        if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
            if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
                IDarpos[0]=0;
                for (int i=1; i<nploc; i++) {
                    IDarpos[i] =IDarpos[i-1]+nIDar_each[i-1];
                }
            }
        }
    }
    
    
    IDuns = new int[Ng];    // array storing the IDs of delete events in current fix (hence local, in the sense that it is both local of the subcomm but also local to current fix_delete)
    if (key>0) {
        int dest = 0;
        MPI_Send(&tIDarr[0], naP, MPI_INT, dest, 1, (universe->subcomm));
    }
    if (key==0) {
        for (int i=0; i<nID_each[0]; i++) {
            IDuns[i] = tIDarr[i];
        }
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&IDuns[IDpos[source]], nID_each[source], MPI_INT, source, 1, (universe->subcomm), &status);
        }
        
        
        // print assembled IDuns on submaster
        if (msk->wplog) {
            std::string msg = "\nNumber of atoms in this fix: ";
            std::ostringstream ss;    ss << Ng;   msg = msg+ss.str();
            output->toplog(msg);
            msg="";   ss.str("");   ss.clear();
            output->toplog("\npos IDuns");
            for (int i=0; i<Ng; i++){
                ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
                ss << IDuns[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                output->toplog(msg);
                msg="";
            }
        }
    }
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            
            if (key>0) {
                int dest = 0;
                MPI_Send(&IDar[0], nR, MPI_INT, dest, 1, (universe->subcomm));
                MPI_Send(&Rar[0], nR, MPI_DOUBLE, dest, 1, (universe->subcomm));
            }
            if (key==0) {
                
                totRp = 0;  // total number of real particles across all procs in this fix
                for (int i=0; i<nploc; i++) totRp += nIDar_each[i];
                
                IDaruns = new int[totRp];    // array storing the IDs of rael particles in current fix (hence local, in the sense that it is both local of the subcomm but also local to current fix_delete)
                Raruns = new double[totRp];    // array storing the IDs of rael particles in current fix (hence local, in the sense that it is both local of the subcomm but also local to current fix_delete)
                
                for (int i=0; i<nIDar_each[0]; i++) {
                    IDaruns[i] = IDar[i];
                    Raruns[i] = Rar[i];
                }
                for (int source=1; source<nploc; source++) {
                    MPI_Recv(&IDaruns[IDarpos[source]], nIDar_each[source], MPI_INT, source, 1, (universe->subcomm), &status);
                    MPI_Recv(&Raruns[IDarpos[source]], nIDar_each[source], MPI_DOUBLE, source, 1, (universe->subcomm), &status);
                }
                
                
                // print assembled ID and Radii arrays of all real particles on submaster
                if (msk->wplog) {
                    std::string msg = "\nNumber of real atoms in this fix: ";
                    std::ostringstream ss;    ss << totRp;   msg = msg+ss.str();
                    output->toplog(msg);
                    msg="";   ss.str("");   ss.clear();
                    output->toplog("\npos IDaruns Raruns");
                    for (int i=0; i<totRp; i++){
                        ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
                        ss << IDaruns[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        ss << Raruns[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                        output->toplog(msg);
                        msg="";
                    }
                }
            }
        }
    }
}



// ---------------------------------------------------------------
// assembles all pair arrays from each processor into block-wise array in submaster, where each block is the locLMP and aDIST in a processor
void Fix_delete::pair_arr_to_submaster(int pos)
{
    
   
    //MPI send size of each local array to submaster. Submaster creates array with enough space and assigns positions to accept local arrays.  All procs then send local arrays to submaster
    nlocR_each = new int[nploc];
    nlocR_each[key] = nlocR;
    
    if (key>0) {
        int dest = 0;
        MPI_Send(&nlocR_each[key], 1, MPI_INT, dest, 1, (universe->subcomm));
    }
    if (key==0 && nploc>1) {
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&nlocR_each[source], 1, MPI_INT, source, 1, (universe->subcomm), &status);
        }
    }
    
    
    // pass local arrays to submaster
    SARpos = new int[nploc];  // position of local array in submaster's  array
    if (key==0) {
        SARpos[0]=0;
        for (int i=1; i<nploc; i++) {
            SARpos[i] =SARpos[i-1]+nlocR_each[i-1];
        }
    
        // allocating array storing the local arrays for the interactions in current fix
        int allrows = SARpos[nploc-1]+nlocR_each[nploc-1];
        int nbytes = ((int) sizeof(double)) * 4 * allrows;
        double *data = (double *) malloc(nbytes);
        nbytes = ((int) sizeof(double *)) * allrows;
        SAR = (double **) malloc(nbytes);
        
        int n = 0;
        for (int i = 0; i < allrows ; i++) {
            SAR[i] = &data[n];
            n += 4;
        }
        Dsub = new double[allrows];
    }
    
    
    if (key>0) {
        int dest = 0;
        MPI_Send(&locLMP[0][0], 4.*nlocR , MPI_DOUBLE, dest, 3, (universe->subcomm));
        MPI_Send(&aDIST[0], nlocR , MPI_DOUBLE, dest, 4, (universe->subcomm));
    }
    if (key==0) {
        for (int j=0; j<nlocR; j++){
            for (int i=0; i<4; i++) {
                SAR[j][i] = locLMP[j][i];
            }
            Dsub[j] = aDIST[j];
        }
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&SAR[SARpos[source]][0], 4.*nlocR_each[source], MPI_DOUBLE, source, 3, (universe->subcomm), &status);
            MPI_Recv(&Dsub[SARpos[source]], nlocR_each[source], MPI_DOUBLE, source, 4, (universe->subcomm), &status);
        }
        
        // print assembled SAR on submaster
        if (msk->wplog) {
            std::string msg = "\nNumber of local pairs: ";
            std::ostringstream ss;    ss << SARpos[nploc-1]+nlocR_each[nploc-1];   msg = msg+ss.str();
            output->toplog(msg);
            msg="";   ss.str("");   ss.clear();
            output->toplog("\npos lmp_pos id1 id2 type1 type2 dist");
            //sleep(1);
            for (int i=0; i<SARpos[nploc-1]+nlocR_each[nploc-1]; i++){
                //sleep(1);
                ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
                ss << SAR[i][0];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                ss << SAR[i][1];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                ss << SAR[i][2];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                ss << SAR[i][3];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                ss << Dsub[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                output->toplog(msg);
                msg="";
            }
        }
    }
    

}




// ---------------------------------------------------------------
// submaster sorts the IDs of all particles in all processors for current fix, and puts them at the end of the IDsrt vector that is for all delete fixes in the subcomm
void Fix_delete::submaster_sort_IDs(int pos)
{
    int posit;

    // the first ID is just added to the vector
    if (Ng>0) {
        IDsrt.push_back(IDuns[0]);
        rates.push_back(0.);
    }
    
    //StoU.push_back(0);  // just create map vectors to be used in submaster_sort_Ids function
    //UtoS.push_back(0);
    // for all subsequent IDs
    for (int j=1; j<Ng; j++) {
        //StoU.push_back(0); // just create map vectors to be used in submaster_sort_Ids function
        //UtoS.push_back(0);
        // if ID is the largest, put it at the end
        if (IDuns[j]>IDsrt[IDsrt.size()-1]) IDsrt.push_back(IDuns[j]);
        else {
            //otherwise look for a position betwwen a smaller and large one
            if (IDuns[j]<IDsrt[fix->fKMCfirst[pos]]) posit = fix->fKMCfirst[pos];
            else {
                //binary search
                int pre = fix->fKMCfirst[pos];
                int post=IDsrt.size()-1;
                while (pre < post) {
                    posit = (int)((pre+post)/2.);
                    if (IDsrt[posit] < IDuns[j]) pre=posit+1;
                    else post=posit;
                }
                posit=pre;
            }
            // shift IDsrt
            IDsrt.push_back(IDsrt[IDsrt.size()-1]);   //shif the last element one up
            for (int k=IDsrt.size()-1; k>posit; k--) {
                IDsrt[k] = IDsrt[k-1];  // shift all other elements until posit
            }
            IDsrt[posit]=IDuns[j]; // record aID in posit
        }
        rates.push_back(0.);  // here is where the rates container is initialised, every time a new IDsrt is added
    }
    
    
    // sorting IDs and corresponding radii for all real particles if a micro mechanism is used
    // the first ID is just added to the vector
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            
            // IDs and radii of all real particles are only used to compute coverage fractions. No need to remember them from one fix delete to the next one (unlike IDsrt corresponding to rates above)
            IDarsrt.clear();
            Rarsrt.clear();
            
            // write down first ID and R
            if (totRp>0) {
                IDarsrt.push_back(IDaruns[0]);
                Rarsrt.push_back(Raruns[0]);
            }
            
            // for all subsequent IDs
            for (int j=1; j<totRp; j++) {
                // if ID is the largest, put it at the end
                if (IDaruns[j]>IDarsrt[IDarsrt.size()-1]) {
                    IDarsrt.push_back(IDaruns[j]);
                    Rarsrt.push_back(Raruns[j]);
                }
                else {
                    //otherwise look for a position betwwen a smaller and large one
                    if (IDaruns[j]<IDarsrt[0]) posit = 0;
                    else {
                        //binary search
                        int pre = 0;
                        int post=IDarsrt.size()-1;
                        while (pre < post) {
                            posit = (int)((pre+post)/2.);
                            if (IDarsrt[posit] < IDaruns[j]) pre=posit+1;
                            else post=posit;
                        }
                        posit=pre;
                    }
                    // shift IDarsrt and Rarsrt
                    IDarsrt.push_back(IDarsrt[IDarsrt.size()-1]);   //shif the last element one up
                    Rarsrt.push_back(Rarsrt[Rarsrt.size()-1]);
                    for (int k=IDarsrt.size()-1; k>posit; k--) {
                        IDarsrt[k] = IDarsrt[k-1];  // shift all other elements until posit
                        Rarsrt[k] = Rarsrt[k-1];
                    }
                    IDarsrt[posit]=IDaruns[j]; // record aID in posit
                    Rarsrt[posit]=Raruns[j];
                }
            }
        }
    }
    
    if (msk->wplog) {
        std::string msg;
        std::ostringstream ss;
        ss << "Total number of real particles is: " << totRp << "\n";
        msg = ss.str(); ss.str("");   ss.clear();
        output->toplog(msg);
        msg = "pos IDarsrt Rarsrt";
        output->toplog(msg);
        for (int i=0; i<totRp; i++) {
            ss << i<<" "<<IDarsrt[i] <<" "<<Rarsrt[i]<<"";
            msg = ss.str(); ss.str("");   ss.clear();
            output->toplog(msg);
        }
    }
    
    
    
}






// ---------------------------------------------------------------
//  map IDs positions from IDuns to IDsrt (the former with size Ng, the latter with size comprising all particles for all delete fixes on that subcom)
void Fix_delete::submaster_map_ID(int pos)
{
    // find ID position in existing sorted vector (portion of IDsrt corresponding to current event, knowing first position id from fix-fMKCfirst and number of events Ng)
    int posit;
    UtoS.clear();
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            StoU.clear();
            StoU.resize(Ng);   // set map vector size
        }
    }
    
    for (int j=0; j<Ng; j++) {
        //binary search
        int pre = fix->fKMCfirst[pos];
        int post=fix->fKMCfirst[pos]+Ng-1;
        while (pre < post) {
            posit = (int)((pre+post)/2.);
            if (IDsrt[posit] < IDuns[j]) pre=posit+1;
                else post=posit;
            }
        posit=pre;
        UtoS.push_back(posit);
        if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
            if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
                StoU[posit-fix->fKMCfirst[pos]] = j;
            }
        }
    }

    
    // submaster prints assembled IDuns, map from unsorted to sorted, and sorted IDs
    if (msk->wplog) {
        std::string msg = "\nNumber of atoms in fix "+fix->fKMCname[pos]+" :";
        std::ostringstream ss;    ss << Ng;   msg = msg+ss.str();
        output->toplog(msg);
        msg="";   ss.str("");   ss.clear();
        output->toplog("\npos IDuns UtoS");
        for (int i=0; i<Ng; i++){
            ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
            ss << IDuns[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << UtoS[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            output->toplog(msg);
            msg="";
        }
        output->toplog("\nWhole sorted ID vector up to current fix \npos IDsrt");
        for (int i=0; i<IDsrt.size(); i++){
            ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
            ss << IDsrt[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            output->toplog(msg);
            msg="";
        }
        if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
            if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
                output->toplog("\nFix-specific sorted-to-unsorted map\npos IDsrt StoU");
                for (int i=0; i<Ng; i++){
                    ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
                    ss << IDsrt[fix->fKMCfirst[pos]+i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                    ss << StoU[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
                    output->toplog(msg);
                    msg="";
                }
            }
        }
    }
}




// ---------------------------------------------------------------
//  The submaster computes coverage areas of particles in the current fix and communicates them back to the other processors in its subcomm
void Fix_delete::submaster_comp_cover(int pos)
{
    CFuns = new double[Ng];
    GMuns = new double[Ng];
    fGMuns = new bool[Ng];
    
    for (int i=0; i<Ng; i++) {
        CFuns[i]=0.;
        GMuns[i]=0.;
        fGMuns[i]=true;
    }
    
    if (msk->wplog) {
        std::string msg = "Number of pairs in neighbor array for "+fix->fKMCname[pos]+" :";
        std::ostringstream ss;    ss << SARpos[nploc-1]+nlocR_each[nploc-1];   msg = msg+ss.str();
        output->toplog(msg);
        msg="";   ss.str("");   ss.clear();
        // pring ids and positions in unsorted atom arrays
        output->toplog("npos id1 id2 type1 type2 upos1 upos2");
    }
    
    
    // for all the rows in the assembled arrays in the submaster...
    for (int i=0; i<SARpos[nploc-1]+nlocR_each[nploc-1]; i++){
        int id1 = SAR[i][0];
        int id2 = SAR[i][1];
        int t1 = SAR[i][2];
        int t2 = SAR[i][3];
        bool flag_t1 = (t1 == fix->fKMCptype[pos]);
        bool flag_t2 = (t2 == fix->fKMCptype[pos]);
        int up1,up2;   // positions of particles 1 and 2 in pair in unsorted vectors (ID and radii)
        up1 = -1;
        up2 = -1;
        
        if (flag_t1){
            // find first particle position in unsorted IDs vector
            int pre = fix->fKMCfirst[pos];
            int post = fix->fKMCfirst[pos]+Ng-1;
            int posit;
            while (pre < post) {
                posit = (int)((pre+post)/2.);
                if (IDsrt[posit] < id1) pre=posit+1;
                    else post=posit;
                }
            posit=pre;
            up1 = StoU[posit - fix->fKMCfirst[pos]];
        }
        
        if (flag_t2){
            // find second particle position in unsorted IDs vector
            int pre = fix->fKMCfirst[pos];
            int post = fix->fKMCfirst[pos]+Ng-1;
            int posit;
            while (pre < post) {
                posit = (int)((pre+post)/2.);
                if (IDsrt[posit] < id2) pre=posit+1;
                    else post=posit;
                }
            posit=pre;
            up2 = StoU[posit - fix->fKMCfirst[pos]];
        }
        
        if (msk->wplog) {
            std::string msg;
            std::ostringstream ss;
            ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
            ss << id1;   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << id2;   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << SAR[i][2];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << SAR[i][3];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << up1;   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << up2;   msg += ss.str()+" ";  ss.str("");   ss.clear();
            output->toplog(msg);
        }
        
       
        if (flag_t1 || flag_t2){
            // compute coverage area, used later to compute number of layers to dissolve in micro particle
            // coverage is given by contact cross section weighted by a distance-dependent factor; the latter is user-provided through two thresholds, e0 below which contact is 100%, ef above which contact is 0%. Linear interpolation in between
            
            double Aij;    // cross section of inter-particle contact
            double Ri,Rj;   // radii of interaction particles in current pair
            double Rij;     // harmonic average of radii of particles in contact
            int rp1,rp2;    // position of interacting particles in the "all real" sorted vectors
            
            // find first particle position in sorted IDar vector
            int pre = 0;
            int post = IDarsrt.size()-1;
            int posit;
            while (pre < post) {
                posit = (int)((pre+post)/2.);
                if (IDarsrt[posit] < id1) pre=posit+1;
                    else post=posit;
                }
            posit=pre;
            rp1 = posit;
            
            //fprintf(screen,"DEBUG 0: rp1 = %d\n",rp1);
            
            
            
            // find second particle position in sorted IDar vector
            pre = 0;
            post = IDarsrt.size()-1;
            while (pre < post) {
                posit = (int)((pre+post)/2.);
                if (IDarsrt[posit] < id2) pre=posit+1;
                    else post=posit;
                }
            posit=pre;
            rp2 = posit;
            

            
            Ri = Rarsrt[rp1];
            Rj = Rarsrt[rp2];
            

            
            if (msk->wplog) {
                std::string msg;
                std::ostringstream ss;
                ss << "rp1 = " << rp1 << "; rp2 = "<< rp2 << "; Ri = " << Ri << "; Rj = " << Rj;
                msg = ss.str();  ss.str("");   ss.clear();
                output->toplog(msg);
            }
            
            Rij = 2. * Ri * Rj/( Ri + Rj);
            Aij = M_PI * Rij * Rij;
            
            // weigh the contact cross section by the distance
            double Dij; // arithmetic average of diameters in contact
            Dij = Ri + Rj;
            double efij = chem->ef[t1-1][t2-1];
            double e0ij = chem->e0[t1-1][t2-1];
            
            if (Dsub[i] > efij * Dij)        Aij = 0.;
            else if (Dsub[i] > e0ij * Dij)   Aij *= (efij - Dsub[i]/Dij) / (efij- e0ij);

            // if first atom is of correct type for this fix, add contact area to the contact fraction arrays (still areas here; will be converted to area fractions later on below)
            if (flag_t1)  {
                CFuns[up1] += Aij/(M_PI*4.*Ri*Ri);
                if (Aij > 0.){
                    double tempGM = 0.;
                    if (t1 != t2) tempGM = - chem->gij[t1-1][t1-1] - chem->gij[t1-1][t2-1] + chem->gij[t2-1][t2-1];
                    if (fGMuns[up1]) {
                        GMuns[up1] = tempGM;
                        fGMuns[up1] = false;
                    }
                    else if (tempGM > GMuns[up1]) GMuns[up1] = tempGM;
                }
            }
            
            if (flag_t2){
                CFuns[up2] += Aij/(M_PI*4.*Rj*Rj);;
                if (Aij > 0.){
                    double tempGM = 0.;
                    if (t1 != t2) tempGM = - chem->gij[t2-1][t2-1] - chem->gij[t2-1][t1-1] + chem->gij[t1-1][t1-1];
                    if (fGMuns[up2]) {
                        GMuns[up2] = tempGM;
                        fGMuns[up2] = false;
                    }
                    else if (tempGM > GMuns[up2]) GMuns[up2] = tempGM;
                }
            }
        }
        
    }
    double nbulk = 12.;   // number of particles in the bulk for a monodisperse system of particls with type selected for dissolution. In future implementations, this could be passed as an input when calling the fix_delete (becasue it applies to all particles to be considered for deletion here)
    for (int i=0; i<Ng; i++){
        CFuns[i] *= 4./nbulk;    // this is a correction such that a speherical particle in contact with 12 equally sized neighbours in an FCC lattice ends up having coverage fraction = 1   (If we did not divide by 3, such a bulk particle would have a coverage fraction of 300%)
    }
    
    
    if (msk->wplog) {
        std::string msg = "\nPer-particle coverage fraction";
        output->toplog(msg);
        msg="";
        // pring ids and covrage fractions of unsorted atoms
        output->toplog("npos id CFuns GMuns");
        std::ostringstream ss;
        for (int i=0; i<Ng; i++){
            ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
            ss << IDuns[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << CFuns[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << GMuns[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            output->toplog(msg);
        }
        
    }
    
}



// ---------------------------------------------------------------
//  each processor in the subcommunicator receives its chunk of unsorted coverage area fractions from the submaster
void Fix_delete::cover_from_submaster(int pos)
{
    
    if (key==0) {
        for (int dest=1; dest<nploc; dest++) {
            MPI_Send(&CFuns[IDpos[dest]], nID_each[dest], MPI_DOUBLE, dest, 1, (universe->subcomm));
            MPI_Send(&GMuns[IDpos[dest]], nID_each[dest], MPI_DOUBLE, dest, 1, (universe->subcomm));
        }
        tCF = new double[nID_each[0]];
        tGM = new double[nID_each[0]];
        for (int i=0; i<nID_each[0]; i++){
            tCF[i] = CFuns[i];
            tGM[i] = GMuns[i];
        }
    }
    else{
        tCF = new double[nID_each[key]];
        tGM = new double[nID_each[key]];
        int source = 0;
        MPI_Recv(&tCF[0], nID_each[key], MPI_DOUBLE, source, 1, (universe->subcomm), &status);
        MPI_Recv(&tGM[0], nID_each[key], MPI_DOUBLE, source, 1, (universe->subcomm), &status);
    }
    
    // Each processor prints its coverage fraction array to its log file
    if (msk->wplog) {
        std::string msg = "\nNumber of local coverage area fractions: ";
        std::ostringstream ss;    ss << nID_each[key];   msg = msg+ss.str();
        output->toplog(msg);
        msg="";   ss.str("");   ss.clear();
        output->toplog("\npos id tCF tGM");
        //sleep(1);
        for (int i=0; i<nID_each[key]; i++){
            //sleep(1);
            ss << i;        msg = ss.str()+" ";   ss.str("");   ss.clear();
            ss << tID[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << tCF[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            ss << tGM[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            output->toplog(msg);
            msg="";
        }
    }
    
}






// ---------------------------------------------------------------
//  each processor computes its deletion rates
void Fix_delete::comp_rates_allpar(int pos)
{
    
    // for all particles in current ("each") processor
    for (int i =0; i<nID_each[key]; i++) {
        
        if (msk->wplog) {
            std::string msg = "\nPARTICLE ";
            std::ostringstream ss;    ss << tID[i];     msg = msg+ss.str()+" ";
            ss.str("");     ss.clear();   ss << tR[i];  msg = msg+ss.str();
            output->toplog(msg);
        }
        
       
        double Pv = 4./3. * M_PI * tR[i] * tR[i] * tR[i];  // particle volume  FOR SPHERICAL
        double Ps = 4. * M_PI * tR[i] * tR[i];  // particle surface   FOR SPHERICAL
        
        int nrt = 1;    //number of reaction in series in chain
        double nrv = 1;    // number of unit reactions or unit chains in particle volume
        
        int chID = -1;  // chain ID: will be > -1 only if mechanism calls indeed a chain
        int rxid = -1;  // reaction ID: will be used now but later in loop too
        
        bool errflag = false;   // error to spot if wrong volume changes defined for reactions: see below
        double fgdV = 0;
        
        if (chem->mechchain[mid]) {  //if reaction is a chain
            chID = chem -> mechrcID[mid];
            nrt = (chem->ch_rxID[chID]).size();
            if (chem -> ch_dVp_fgd[chID] < 0) {
                // default for allpar mechanism is linear changes in surface and energy with number of molecules in particle volume, i.e. int_lin
                nrv = round( - Pv / chem -> ch_dVp_fgd[chID]);   // minus because a proper dissolution reax should have negative fgd changes in chemDB
            } else {errflag = true; fgdV=chem -> ch_dVp_fgd[chID];}
        }
        else{ //if instead it is a single reaction
            nrt = 1;
            rxid = chem->mechrcID[mid];
            if (chem -> rx_dVp_fgd[rxid] < 0) {
                nrv = round( - Pv / chem -> rx_dVp_fgd[rxid]);
            } else {errflag = true; fgdV=chem -> rx_dVp_fgd[rxid];}
        }
        
        
        if (errflag == true) {
            std::ostringstream wd1,wd2,wd3,wd4;
            wd1<<me; wd2<<universe->color; wd3<<msk->tempo,wd4<<fgdV;
            std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: fix "+fix->fKMCname[pos]+" is of deletion type, which requires negative fgd volume changes, whereas its associated mechanism has fgd volume change = "+wd4.str()+" \n****************************************************\n***********************************************\n";
            fprintf(screen,"%s",msg.c_str());
        }
        
        std::string msg;
  
        
        
        
        
         // change of surface energy and interaction energy per reaction unit (single reaction or chain: to be further subdivided per step in chain later on)
        double DSpu;  // negative becasue delation removes energy
        double DUpu;;  // change in energy caused by removing the particle
        
        if (strcmp(chem->mechinter[mid].c_str(),"int_1lin")==0)  {
            DSpu = -Ps/pow((double)nrv,1./3.);
            DUpu = -2.*tE[i]/pow((double)nrv,1./3.);
        }
        else if (strcmp(chem->mechinter[mid].c_str(),"int_2lin")==0) {
            DSpu = -Ps/pow((double)nrv,68.5/100.);      // power artificailly tuned to get hetero nucleation in cg sim... just to fiddle around: do not use those in a proper simulation
            DUpu = -2.*tE[i]/pow((double)nrv,68.5/100.);
        }
        else{
            DSpu = -Ps/((double)nrv);
            DUpu = -2.*tE[i]/((double)nrv);
        }
        
        
        // number of repetitions of reaction of chain to delete particle. For "allpar" just 1.
        int nrx = 1;
        
        // compute rate looking at each reaction in the mechanism sequence
        double ri=0. , DTi = 0., DTtot = 0.;  // rate of each reaction in sequence, associated time increement and total cumulative time increment
        
        
        for (int j=0; j<nrx; j++) {
            
            if (msk->wplog) { msg += "; rx "; std::ostringstream ss;    ss << j;   msg = msg+ss.str(); }
            
            // compute rate of reaction sequence
            for (int k=0; k< nrt; k++) { //all the reaction in series in chain seq.
                
                if (msk->wplog) { msg += "; rx_step ";    std::ostringstream ss;    ss << k;   msg += ss.str(); }
                
                // find reaction id, either from chain or the single one
                if (chem->mechchain[mid]) rxid = chem->ch_rxID[chID][k];
                else rxid = chem->mechrcID[mid];   // probably redundant as already defined above: to be checked
                
                double DGx = chem -> compDGx(rxid); // reaction specific
                double cx = chem -> cx[chem->rx_DGID[rxid]];
                double dim = chem -> dim[chem->rx_DGID[rxid]];
                double gammax = chem -> compgammax(rxid);
                double KT = msk->kB * solution->Temp;
                
                if (msk->wplog) {
                    std::ostringstream ss;
                    msg += "; DG* ";    ss << DGx;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; gamma* ";    ss << gammax;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; KB ";    ss << (msk->kB);   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; T ";    ss << (solution->Temp);   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                double kappa = 1.;       // transmission coefficient
                double Sen = chem->compSen(mid);
                
            
                //bool net = false; // the beta calculator in simulation.cpp uses stoichio coefficients that are defined by used in chemDB for forward reaction. If beta is computed for the net rate, the reaction is backward hence stoichio coeff need multiplication by -1
                //double beta = 1.; // beta in prefactor: for "allpar" one should take 1 as per TST. Maybe in future a calculator for this too.
                std::string topass = "reac";
                double Qreac = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                
                double Vt = - (chem -> rx_dVt_fgd[rxid]);  // Tributary volume of fgd deleted by the reaction at the current step in the chain. Minus sign because dV of fdg < 0 for a proper dissolution reaction.
                
               //Delta surface and Delta U are assumed to be proportional to the number of reaction/chain units and relative importance of each step in chain
                double DSi = DSpu;
                if (chem->mechchain[mid]) DSi *= (chem->ch_rdVp_fgd[chID][k]);
                
                double DUi = 0.;
                if (strcmp(chem->mechinter[mid].c_str(),"int_no")!=0) {
                    DUi = DUpu;
                    if (chem->mechchain[mid]) DUi *= (chem->ch_rdV_fgd[chID][k]);
                }
                
                double r0 = kappa * KT / msk->hpl / gammax * cx * exp(-DGx / KT);
                
                if (msk->wplog) {
                    std::ostringstream ss;
                    msg += "; r0 ";    ss << r0;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Qreac ";    ss << Qreac;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Sen ";    ss << Sen;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; DSi ";    ss << DSi;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; DUi ";    ss << DUi;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Vti ";    ss << Vt;   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                
               // double fa = 1.;
                //if (dim==2 || dim==1) fa=3.;
                
                r0 = r0*pow(Vt,dim/3.);
                
                double ki = (chem -> ki[rxid]);
                
                // Rate equation as per TST (see notes_on_TST.pdf by Masoero, 2019, where ki=0 was assumed as per classical TST)
                ri = r0 * Qreac * exp((1.-ki)*(-Sen * DSi - DUi)/KT ) ;
                
                
                // if net rate is requested by user in chemDB, then add rate backward
                if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                    
                    topass = "prod";
                    double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    ri -= r0 * Qprod / chem->Keq[rxid]     * exp(ki*(Sen * DSi + DUi)/KT );
                    
                    if (msk->wplog) {
                        std::ostringstream ss;
                        msg += "; ri ";    ss << ri;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Qprod ";    ss << Qprod;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Keq ";    ss << chem->Keq[rxid] ;   msg += ss.str();    ss.str(""); ss.clear();
                    }
                    
                }
                
                
                
                if (ri < 0.) ri = 0.;
                
                DTi = 1./ri;
                DTtot += DTi;
                
                if (msk->wplog) {
                    msg += ", DT ";
                    std::ostringstream ss;    ss << DTi;   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", DTtot ";
                    ss << DTtot;   msg += ss.str();
                }
            }
            
            if (msk->wplog) {
                msg += "\n";
                output->toplog(msg);
            }
        }
        rate_each[i] = 1./DTtot;
        if (rate_each[i]<0.) rate_each[i] = 0.;    // if the backward ri's are > forward ri's the overall rate may end up < 0, which means that the current deletion event should not happen, hence its rate should be zero (not negative..)   CHECK THAT THIS DOES NOT GIVE PROBLEMS WHEN SELECTING THE EVENT TO CARRY OUT FROM CUMULATIVE RATE VECTORS, IN PARTICULAR WITH THE BINARY SEARCH ALGORITHM
    
        
        if (msk->wplog) {
            std::ostringstream ss;
            std::string msg = "\nPARTICLE ";
            ss << tID[i];   msg += ss.str();  ss.str("");   ss.clear();
            msg += " has deletion rate  ";
            ss << rate_each[i];   msg += ss.str();
            output->toplog(msg);
        }
    }
    
}




// ---------------------------------------------------------------
//  each processor computes its deletion rates
void Fix_delete::comp_rates_micro(int pos)
{
    
    // for all particles in current processor, compute the dissolution rate
    for (int i=0; i<nID_each[key]; i++) {
        
        if (msk->wplog) {
            std::string msg = "\n -----------------------------------\nPARTICLE ";
            std::ostringstream ss;    ss << tID[i];     msg = msg+ss.str()+" ";
            ss.str("");     ss.clear();   ss << tR[i];  msg = msg+ss.str()+" ";
            ss.str("");     ss.clear();   ss << tCF[i];  msg = msg+ss.str()+" ";
            ss.str("");     ss.clear();   ss << tGM[i];  msg = msg+ss.str()+"\n";
            output->toplog(msg);
        }
        
        double Pv = 4./3. * M_PI * tR[i] * tR[i] * tR[i];  // particle volume  FOR SPHERICAL
        //double Ps = 4. * M_PI * tR[i] * tR[i];  // particle surface   FOR SPHERICAL - this was used to compute changes in surface area per molecule, but in this micro mechanism, surface energy is not considered (particle dissolution is assumed to be kink-driven)
        
        int nrt = 1;    //number of reaction in series in chain
        double nrv = 1;    // number of unit reactions or unit chains in particle volume
        
        int chID = -1;  // chain ID: will be > -1 only if mechanism calls indeed a chain
        int rxid = -1;  // reaction ID: will be used now but later in loop too
        
        bool errflag = false;   // error to spot if wrong volume changes defined for reactions: see below
        double fgdV = 0;
        
        if (chem->mechchain[mid]) {  //if reaction is a chain
            chID = chem -> mechrcID[mid];
            nrt = (chem->ch_rxID[chID]).size();
            if (chem -> ch_dVp_fgd[chID] < 0) {
                // default for allpar mechanism is linear changes in surface and energy with number of molecules in particle volume, i.e. int_lin
                nrv = round( - Pv / chem -> ch_dVp_fgd[chID]);   // minus because a proper dissolution reax should have negative fgd changes in chemDB
            } else {errflag = true; fgdV=chem -> ch_dVp_fgd[chID];}
        }
        else{ //if instead it is a single reaction
            nrt = 1;
            rxid = chem->mechrcID[mid];
            if (chem -> rx_dVp_fgd[rxid] < 0) {
                nrv = round( - Pv / chem -> rx_dVp_fgd[rxid]);
            } else {errflag = true; fgdV=chem -> rx_dVp_fgd[rxid];}
        }
        
        
        if (errflag == true) {
            std::ostringstream wd1,wd2,wd3,wd4;
            wd1<<me; wd2<<universe->color; wd3<<msk->tempo,wd4<<fgdV;
            std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: fix "+fix->fKMCname[pos]+" is of deletion type, which requires negative fgd volume changes, whereas its associated mechanism has fgd volume change = "+wd4.str()+" \n****************************************************\n***********************************************\n";
            fprintf(screen,"%s",msg.c_str());
        }
  
        
        
         // change of  interaction energy per reaction unit (single reaction or chain: to be further subdivided per step in chain later on)
        //double DSpu;  // negative becasue delation removes energy
        double DUpu;  // change in energy caused by removing the particle
        
        if (strcmp(chem->mechinter[mid].c_str(),"int_1lin")==0)  {
            //DSpu = -Ps/pow((double)nrv,1./3.);
            //DUpu = -2.*tE[i]/pow((double)nrv,1./3.);
        }
        else if (strcmp(chem->mechinter[mid].c_str(),"int_2lin")==0) {
            //DSpu = -Ps/pow((double)nrv,68.5/100.);      // power artificailly tuned to get hetero nucleation in cg sim... just to fiddle around: do not use those in a proper simulation
            //DUpu = -2.*tE[i]/pow((double)nrv,68.5/100.);
        }
        else{
            //DSpu = -Ps/((double)nrv);
            DUpu = -2.*tE[i]/((double)nrv);
        }
        

        // number of layers to dissolve (it will depend on fractional coverage area)
        double lim_CF= std::stod(chem->mechpar[mid][3]); //limiting coverage fraction
        //fprintf(screen,"\n The value of limiting coverage fraction is %f\n",lim_CF);

        int nrL = 0;
        
        double unit_thick;    // thickness of dissolving unit in radial direction of particle
        if (chem->mechchain[mid]) {  //if reaction is a chain
            unit_thick = pow(- chem->ch_dVp_fgd[chID],1./3.);
        }
        else{ //if instead it is a single reaction
            unit_thick = pow(- chem->rx_dVp_fgd[rxid],1./3.);
        }
        
        if (tCF[i] < 0.5) nrL = round(tR[i]/unit_thick);
        else if (tCF[i] <= lim_CF) nrL = round(2.*tR[i]/unit_thick);
        
        bool flag_bulk = false;
        if (nrL == 0) flag_bulk = true;
        
        
        
        
        // number of units to dissolve a full layer (depends on density of kinks per layer)
        int nrS = 1;
        if (chem->mechchain[mid])  nrS = round(1./chem->ch_Fk[chID]);
        else nrS = round(1./chem->rx_Fk[rxid]);
           

        
        
        
        
        // compute rates for all reaction units (multi-step if a chain), for all layers, for all surface stpes
        
        double ri=0. , DTi = 0., DTtot = 0.;  // rate of each reaction in sequence, associated time increement and total cumulative time increment
        
        for (int j=0; j<nrL-1; j++){

            std::string msg;
            if (msk->wplog) {
                msg += "PARTICLE ";
                std::ostringstream ss;    ss << i;   msg = msg+ss.str()+"; Layer ";
                ss.str("");     ss.clear();   ss << j+1;  msg = msg+ss.str()+" out of ";
                ss.str("");     ss.clear();   ss << nrL;  msg = msg+ss.str()+"; ";
            }
            
            for (int jj=0; jj<nrS; jj++){
                
                if (msk->wplog) {
                    msg += "On-surface unit number ";
                    std::ostringstream ss;    ss << jj+1;   msg = msg+ss.str()+" out of ";
                    ss.str("");     ss.clear();   ss << nrS;  msg = msg+ss.str()+"; ";
                }
                
                // compute rate of reaction sequence
                for (int k=0; k< nrt; k++) { //all the reaction in series in chain seq.
                    
                    if (msk->wplog) { msg += "; rx_step ";    std::ostringstream ss;    ss << k;   msg += ss.str();}
                    
                    // find reaction id if it is a chain, or keep previous one if a single reax
                    if (chem->mechchain[mid]) rxid = chem->ch_rxID[chID][k];
                    
                    double DGx = chem -> compDGx(rxid); // reaction specific
                    double cx = chem -> cx[chem->rx_DGID[rxid]];
                    double dim = chem -> dim[chem->rx_DGID[rxid]];
                    double gammax = chem -> compgammax(rxid);
                    double KT = msk->kB * solution->Temp;
                    
                    if (msk->wplog) {
                        std::ostringstream ss;
                        msg += "; DG* ";    ss << DGx;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; c* ";    ss << cx;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; gamma* ";    ss << gammax;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; KB ";    ss << (msk->kB);   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; T ";    ss << (solution->Temp);   msg += ss.str();    ss.str(""); ss.clear();
                    }
                    
                    double kappa = 1.;       // transmission coefficient
                    
                    std::string topass = "reac";
                    double Qreac = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    
                    double Vt = - (chem -> rx_dVt_fgd[rxid]);  // Tributary volume of fgd deleted by the reaction at the current step in the chain. Minus sign because dV of fdg < 0 for a proper dissolution reaction.
                    
                    double DUi = 0.;
                    if (strcmp(chem->mechinter[mid].c_str(),"int_no")!=0) {
                        DUi = DUpu;
                        if (chem->mechchain[mid]) DUi *= (chem->ch_rdV_fgd[chID][k]);
                        
                        // Add change in strain energy internal to the molecules
                        for (int j=0; j<chem->fgd_molID[rxid].size(); j++){
                            int molid = chem->fgd_molID[rxid][j];
                            DUi += chem->mol_Us[molid] * chem->mol_vm[molid] * chem->fgd_nmol[rxid][j];
                        }
                    }
                    
                    double r0 = kappa * KT / msk->hpl / gammax * cx * exp(-DGx / KT);
                    
                    if (msk->wplog) {
                        std::ostringstream ss;
                        msg += "; r0 ";    ss << r0;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Qreac ";    ss << Qreac;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; DUi ";    ss << DUi;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Vti ";    ss << Vt;   msg += ss.str();    ss.str(""); ss.clear();
                    }
                    
                    r0 = r0*pow(Vt,dim/3.);
                    
                    double ki = (chem -> ki[rxid]);
                    
                    // Forward rate equation
                    ri = r0 * Qreac * exp((1.-ki)*(- DUi)/KT ) ;
                    
                    
                    // if net rate is requested by user in chemDB, then add rate backward
                    if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                        
                        topass = "prod";
                        double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                        
                        ri -= r0 * Qprod / chem->Keq[rxid] * exp(ki * DUi / KT );
                        
                        if (msk->wplog) {
                            std::ostringstream ss;
                            msg += "; ri ";    ss << ri;   msg += ss.str();    ss.str(""); ss.clear();
                            msg += "; Qprod ";    ss << Qprod;   msg += ss.str();    ss.str(""); ss.clear();
                            msg += "; Keq ";    ss << chem->Keq[rxid] ;   msg += ss.str();    ss.str(""); ss.clear();
                        }
                        
                    }
                    
                    if (ri < 0.) ri = 0.;
                    
                    DTi = 1./ri;
                    DTtot += DTi;
                    
                    if (msk->wplog) {
                        msg += ", DT ";
                        std::ostringstream ss;    ss << DTi;   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", DTtot ";
                        ss << DTtot;   msg += ss.str();
                    }
                }
            }
            if (msk->wplog){
                msg+="\n";
                output->toplog(msg);
            }
        }
        
        
        
        // Add to Dtot the contribution from the last layer
        std::string msg;
        if (msk->wplog) {
            msg += "PARTICLE ";
            std::ostringstream ss;    ss << i;   msg = msg+ss.str()+"; Last layer; ";
        }
        
        for (int jj=0; jj<nrS; jj++){
            
            if (msk->wplog) {
                msg += "On-surface unit number ";
                std::ostringstream ss;    ss << jj+1;   msg = msg+ss.str()+" out of ";
                ss.str("");     ss.clear();   ss << nrS;  msg = msg+ss.str()+"; ";
            }
            
            // compute rate of reaction sequence
            for (int k=0; k< nrt; k++) { //all the reaction in series in chain seq.
                
                if (msk->wplog) { msg += "; rx_step ";    std::ostringstream ss;    ss << k;   msg += ss.str();}
                
                // find reaction id if it is a chain, or keep previous one if a single reax
                if (chem->mechchain[mid]) rxid = chem->ch_rxID[chID][k];
                
                double DGx = chem -> compDGx(rxid); // reaction specific
                double cx = chem -> cx[chem->rx_DGID[rxid]];
                double dim = chem -> dim[chem->rx_DGID[rxid]];
                double gammax = chem -> compgammax(rxid);
                double KT = msk->kB * solution->Temp;
                
                if (msk->wplog) {
                    std::ostringstream ss;
                    msg += "; DG* ";    ss << DGx;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; c* ";    ss << cx;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; gamma* ";    ss << gammax;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; KB ";    ss << (msk->kB);   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; T ";    ss << (solution->Temp);   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                double kappa = 1.;       // transmission coefficient
                
                std::string topass = "reac";
                double Qreac = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                
                double Vt = - (chem -> rx_dVt_fgd[rxid]);  // Tributary volume of fgd deleted by the reaction at the current step in the chain. Minus sign because dV of fdg < 0 for a proper dissolution reaction.
                
                double DUi = 0.;
                if (strcmp(chem->mechinter[mid].c_str(),"int_no")!=0) {
                    DUi = DUpu;
                    if (chem->mechchain[mid]) DUi *= (chem->ch_rdV_fgd[chID][k]);
                    
                    // Add change in strain energy internal to the molecules
                    for (int j=0; j<chem->fgd_molID[rxid].size(); j++){
                        int molid = chem->fgd_molID[rxid][j];
                        DUi += chem->mol_Us[molid] * chem->mol_vm[molid] * chem->fgd_nmol[rxid][j];
                    }
                }
                
                double r0 = kappa * KT / msk->hpl / gammax * cx * exp(-DGx / KT);
                
                if (msk->wplog) {
                    std::ostringstream ss;
                    msg += "; r0 ";    ss << r0;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Qreac ";    ss << Qreac;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; DUi ";    ss << DUi;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Vti ";    ss << Vt;   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                r0 = r0*pow(Vt,dim/3.);
                
                double ki = (chem -> ki[rxid]);
                
                // Forward rate equation: NB this includes contribution from largest interfacial energy at contact with other phases (from function computing coverage areas)
                double uk_area;    // area of a unit kink (assuming cubic units)
                if (chem->mechchain[mid]) {  //if reaction is a chain
                    uk_area = pow(- chem->ch_dVp_fgd[chID],2./3.) * 3.;
                }
                else{ //if instead it is a single reaction
                    uk_area = pow(- chem->rx_dVp_fgd[rxid],2./3.) * 3.;
                }
                
                ri = r0 * Qreac * exp((1.-ki)*( - DUi - tGM[i] * uk_area )/KT ) ;
                
                
                // if net rate is requested by user in chemDB, then add rate backward
                if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                    
                    topass = "prod";
                    double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    ri -= r0 * Qprod / chem->Keq[rxid] * exp(ki * ( DUi + tGM[i] * uk_area)/ KT );
                    
                    if (msk->wplog) {
                        std::ostringstream ss;
                        msg += "; ri ";    ss << ri;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Qprod ";    ss << Qprod;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Keq ";    ss << chem->Keq[rxid] ;   msg += ss.str();    ss.str(""); ss.clear();
                    }
                    
                }
                
                if (ri < 0.) ri = 0.;
                
                DTi = 1./ri;
                DTtot += DTi;
                
                if (msk->wplog) {
                    msg += ", DT ";
                    std::ostringstream ss;    ss << DTi;   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", DTtot ";
                    ss << DTtot;   msg += ss.str();
                }
            }
        }
        if (msk->wplog){
            msg+="\n";
            output->toplog(msg);
        }
        
        
        
        
        if (flag_bulk) rate_each[i] = 0.;
        else rate_each[i] = 1./DTtot;
        
        if (rate_each[i]<0.) rate_each[i] = 0.;    // PROBABLY OBSOLETE ...... if the backward ri's are > forward ri's the overall rate may end up < 0, which means that the current deletion event should not happen, hence its rate should be zero (not negative..)   CHECK THAT THIS DOES NOT GIVE PROBLEMS WHEN SELECTING THE EVENT TO CARRY OUT FROM CUMULATIVE RATE VECTORS, IN PARTICULAR WITH THE BINARY SEARCH ALGORITHM
    
        
        if (msk->wplog) {
            std::ostringstream ss;
            std::string msg = "\nPARTICLE ";
            ss << tID[i];   msg += ss.str();  ss.str("");   ss.clear();
            msg += " has deletion rate  ";
            ss << rate_each[i];   msg += ss.str();
            output->toplog(msg);
        }
        
    }
    
}




// ---------------------------------------------------------------
//  each processor computes its deletion rates
void Fix_delete::comp_rates_allser(int pos)
{
    
    std::ostringstream wd1,wd2,wd3;
    wd1<<me; wd2<<universe->color; wd3<<msk->tempo;
    std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: I have changed fix_delete and not updated the serial mechanism yet, also because I am not sure whether the serial mechanism makes physically sense at all - better mechanisms for coarse graining should be implemented...   If you decide to update the serial mechanism, make sure you pay attention to nrv that, if too big, will overflow the integer. This can happen if using large paritcle, e.g. 1 micrometre in diameter. In such case, the serial mechanism is nonsense. You should therefore put a constraint to the maximum nrv for which it is reasonable to assume a serial mechanism \n****************************************************\n***********************************************\n";
    fprintf(screen,"%s",msg.c_str());
    
   /*
    if (msk->wplog) {
        std::ostringstream ss;
        std::string msg = "\nThis is an all-series mechanism ";
    }
    // for all particles in current ("each") processor
    for (int i =0; i<nID_each[key]; i++) {
        
        if (msk->wplog) {
            std::string msg = "\nPARTICLE ";
            std::ostringstream ss;    ss << tID[i];   msg = msg+ss.str();
            output->toplog(msg);
        }
        
        
        double Pv = 4./3. * M_PI * tR[i] * tR[i] * tR[i];  // particle volume  FOR SPHERICAL
        double Ps = 4. * M_PI * tR[i] * tR[i];  // particle surface   FOR SPHERICAL
        
        int nrt = 1;    //number of reaction in series in chain
        int nrv = 1;    // number of unit reactions or unit chains in particle volume
        
        int chID = -1;  // chain ID: will be > -1 only if mechanism calls indeed a chain
        int rxid = -1;  // reaction ID: will be used now but later in loop too
        if (chem->mechchain[mid]) {  //if reaction is a chain
            chID = chem -> mechrcID[mid];
            nrt = (chem->ch_rxID[chID]).size();
            nrv = round( - Pv / chem -> ch_dV_fgd[chID]);   // minus because a proper dissolution reax should have negative fgd changes defined in the chemDB file
        }
        else{ //if instead it is a single reaction
            nrt = 1;
            rxid = chem->mechrcID[mid];
            nrv = round( - Pv / chem -> rx_dV_fgd[rxid]);
        }
        
        std::string msg;
        if (msk->wplog) { msg = "\n nrv "; std::ostringstream ss;    ss << nrv;   msg = msg+ss.str(); }
        
        // Current particle volume and particle surface to be reduced in following loop every time a reaction occurs
        double Pv_cur = Pv;
        double Ps_cur = Ps;
        
        // number of repetitions of reaction of chain to delete particle. For "allser" it should be equal to number of chains.
        int nrx = nrv;
        
        // compute rate looking at each reaction in the mechanism sequence
        double ri=0. , DTi = 0., DTtot = 0.;  // rate of each reaction in sequence, associated time increment and total cumulative time increment
        
        
        for (int j=0; j<nrx; j++) {
            
            if (msk->wplog) { msg += "; rx "; std::ostringstream ss;    ss << j;   msg = msg+ss.str(); }
            
            // compute rate of reaction sequence
            for (int k=0; k< nrt; k++) { //all the reaction in series in chain seq.
                
                if (msk->wplog) { msg += "; rx_step ";    std::ostringstream ss;    ss << k;   msg += ss.str(); }
                
                // find reaction id, either from chain or the single one
                if (chem->mechchain[mid]) rxid = chem->ch_rxID[chID][k];
                else rxid = chem->mechrcID[mid];   // probably redundant as already defined above: to be checked
                
                double DGx = chem -> compDGx(rxid); // reaction specific
                double gammax = chem -> compgammax(rxid);
                double KT = msk->kB * solution->Temp;
                
                if (msk->wplog) {
                    std::ostringstream ss;
                    msg += "; DG* ";    ss << DGx;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; gamma* ";    ss << gammax;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; KB ";    ss << (msk->kB);   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; T ";    ss << (solution->Temp);   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                double r0star =  KT / msk->hpl / gammax * exp(-DGx / KT);
                double apf = 1.; // molecule size in prefactor: for "allpar" one should take 1 and use gammax above in dimensionless units.
                double ki = chem->ki[rxid];
                double Sen = chem->compSen(mid);
                
                
                bool net = false; // the beta calculator in simulation.cpp uses stoichio coefficients that are defined by used in chemDB for forward reaction. If beta is computed for the net rate, the reaction is backward hence stoichio coeff need multiplication by -1
                double beta = 1.; // beta in prefactor: for "allpar" one should take 1 as per TST. Maybe in future a calculator for this too.
                
                
                // Compute new volume of particle by subtracting 'change in volume due to one reaction' from the current volume; since dissolution causes a negative volume change, a positive sign is used
                double Pv_new = Pv_cur + chem -> rx_dV_fgd[rxid];
                
                // Compute new surface of particle from new particle volume
                double Dnew = pow((Pv_new * 6. / M_PI ), (1./3.));
                
                double Ps_new = (M_PI * pow(Dnew, 2.));
                
                // Compute change in surface due to one reaction
                double DSi = Ps_new - Ps_cur;  // Change in surface due to dissolution is negative, hence current surface is subtracted from new surface
                
                // Compute change in interaction energy due to one reaction
                double DUi = -2.*tE[i] / Ps * (-DSi) ;
                
                // Update particle volume and surface
                Pv_cur = Pv_new ;
                
                Ps_cur = Ps_new ;
                
                
                if (msk->wplog) {
                    std::ostringstream ss;
                    msg += "; r0* ";    ss << r0star;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; apf ";    ss << apf;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; beta ";    ss << beta;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; ki ";    ss << ki;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Sen ";    ss << Sen;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; DSi ";    ss << DSi;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; DUi ";    ss << DUi;   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                
                // Rate equation as per TST (see Shvab et al 2017)
                ri = r0star * apf * apf *  beta * exp( (1. - ki) * (-Sen * DSi - DUi)/KT );
                
                
                // if net rate is requested by user in chemDB, then add rate backward
                if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                    net = true;
                    beta = solution->compbeta(rxid,net,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    ri -= r0star * apf * apf * beta * exp( ki * (Sen * DSi + DUi)/KT );
                }
    
                if (ri < 0.) ri = 0.;
    
                DTi = 1./ri;
                DTtot += DTi;
                
                if (msk->wplog) {
                    msg += ", DT ";
                    std::ostringstream ss;    ss << DTi;   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", DTtot ";
                    ss << DTtot;   msg += ss.str();
                }
            }
            
            if (msk->wplog) {
                msg += "\n";
                output->toplog(msg);
            }
        }
        rate_each[i] = 1./DTtot;
        if (rate_each[i]<0.) rate_each[i] = 0.;    // if the backward ri's are > forward ri's the overall rate may end up < 0, which means that the current deletion event should not happen, hence its rate should be zero (not negative..)   CHECK THAT THIS DOES NOT GIVE PROBLEMS WHEN SELECTING THE EVENT TO CARRY OUT FROM CUMULATIVE RATE VECTORS, IN PARTICULAR WITH THE BINARY SEARCH ALGORITHM
        
        
        if (msk->wplog) {
            std::ostringstream ss;
            std::string msg = "\nPARTICLE ";
            ss << tID[i];   msg += ss.str();  ss.str("");   ss.clear();
            msg += " has rate  ";
            ss << rate_each[i];   msg += ss.str();
            output->toplog(msg);
        }
    }
    */
}





// ---------------------------------------------------------------
//  the submaster receives rates from all slave processors and orders them in the rates vector corresponding to IDsrt (which is the only one that will be remembered)
void Fix_delete::rates_to_submaster(int pos){
   
    unsRates = new double[Ng];    // array storing the rates of delete events in current fix
    if (key>0) {
        int dest = 0;
        MPI_Send(&rate_each[0], naP, MPI_DOUBLE, dest, 1, (universe->subcomm));
    }
    if (key==0) {
        for (int i=0; i<nID_each[0]; i++) {
            unsRates[i] = rate_each[i];
        }
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&unsRates[IDpos[source]], nID_each[source], MPI_DOUBLE, source, 1, (universe->subcomm), &status);
        }
        
        if (msk->wplog) {
            std::string msg = "\nAssembled UNSORTED RATES";
            output->toplog(msg);
            msg = "IDuns   UtoS    unsRate";
            output->toplog(msg);
            for (int i=0; i<Ng; i++) {
                std::ostringstream ss;
                ss << IDuns[i];   msg = ss.str()+"   ";  ss.str("");   ss.clear();
                ss << UtoS[i];   msg += ss.str()+"   ";  ss.str("");   ss.clear();
                ss << unsRates[i];   msg += ss.str();
                output->toplog(msg);
            }
        }
        
        // copy assembled unsorted rate to sorted rate vector containing all fixes of this type on this subcomm
        for (int i=0; i<Ng; i++) {
            rates[UtoS[i]] = unsRates[i];
        }
        
        if (msk->wplog) {
            std::string msg = "\n SORTED RATES:  ";
            std::ostringstream ss;
            ss << rates.size();   msg += ss.str()+" particles";  ss.str("");   ss.clear();
            output->toplog(msg);
            msg = "IDsrt      rates";
            output->toplog(msg);
            for (int i=0; i<rates.size(); i++) {
                ss << IDsrt[i];   msg = ss.str()+"   ";  ss.str("");   ss.clear();
                ss << rates[i];   msg += ss.str(); ss.str("");   ss.clear();
                output->toplog(msg);
            }
        }
        
    }
    delete [] unsRates;
}




// ---------------------------------------------------------------
double Fix_delete::CrateAll()
{

    double Qk = 0.;
    
    for (int i=0; i<rates.size(); i++) Qk += rates[i];
    
    if (msk->wplog) {
        std::string msg = "\n rate vector size is: ";
        std::ostringstream ss;
        ss << rates.size();   msg += ss.str();  ss.str("");   ss.clear();
        output->toplog(msg);
        msg = "\n CUMULATED RATE IS : ";
        ss << Qk;   msg += ss.str();
        output->toplog(msg);
    }
    
    return Qk;
}



// ---------------------------------------------------------------
double Fix_delete::CratesTupdate(double Dt)
{
    double toreturn = 0.;

    if (rates.size() > 0) {
        if (krun->fMC_justreset){
            for (int i=0; i<rates.size(); i++) ratesT.push_back(0.);
            for (int i=0; i<rates.size(); i++) CratesT.push_back(0.);
        }
        
        for (int i=0; i<rates.size(); i++) ratesT[i] += rates[i] * Dt;
        
        CratesT[0] = ratesT[0];
        for (int i=1; i<rates.size(); i++) CratesT[i] = CratesT[i-1] + ratesT[i];
        
        toreturn = CratesT.back();
    }
    
    return toreturn;
}



// ---------------------------------------------------------------
void Fix_delete::findEV(double CRC)
{
    //binary search
    int pre = 0;
    int post = rates.size()-1;
    int posit;
    while (pre < post) {
        posit = (int)((pre+post)/2.);
        if (CratesT[posit] < CRC) pre=posit+1;
        else post=posit;
    }
    posit=pre;
    EVpID = IDsrt[posit];
    
    fprintf(screen,"\n PROC %d posit is %d and IDsrt[posit] is %d:\n",me,posit,IDsrt[posit]);
    //for (int i=0; i<IDsrt.size(); i++) {
      //  fprintf(screen,"\n ID[%d] = %d",i,IDsrt[i]);
    //}
    
    // find absolute fix ID for selected particle
    int i=0;
    bool found_fix = false;
    while (i<cumNPfix.size() && !found_fix){
        if (posit<=cumNPfix[i]-1) {
            found_fix = true;
            EVafixID = fixaID[i];
        }
        else i++;
    }
}




// ---------------------------------------------------------------
// Deleting selected particle with pID
void Fix_delete::execute(int pID)
{
    if (msk->wplog) {
        std::string msg = "\n ****************** \n Deleting particle ";
        std::ostringstream ss;
        ss << pID;   msg += ss.str()+"\n********************\n";
        output->toplog(msg);
    }
    
    
    std::string tolmp = "group sub id ";
    std::ostringstream ss; ss << pID;   tolmp += ss.str();
    lammpsIO->lammpsdo(tolmp);
    
    if (msk->wplog) {
        std::string msg = "\n Succesfully created group in lammps \n";
        output->toplog(msg);
    }
    
    
    // read radius of particle to be deleted
    tolmp = "compute tempRAD sub property/atom radius";
    lammpsIO->lammpsdo(tolmp);
    
    //temp dump to update variables and computes
    tolmp = "dump tdID all custom 1 dump.temp_"+universe->SCnames[universe->color]+" c_tempRAD ";
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "compute craddd sub reduce sum c_tempRAD";
    lammpsIO->lammpsdo(tolmp);
    
    lammpsIO->lammpsdo("run 1");
    
    // alters lammps thermo to compute the required quantities
    tolmp = lammpsIO->lmpThSt + " c_craddd";
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");
    
    tolmp = "variable vidtt equal $(c_craddd)";
    lammpsIO->lammpsdo(tolmp);
    
    // alters lammps thermo to compute the required quantities
    tolmp = lammpsIO->lmpThSt + " c_craddd v_vidtt";
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");
    


    
    
    //tolmp = "compute tempRAD sub property/atom radius";
    //lammpsIO->lammpsdo(tolmp);
    //tolmp = "run 0";
    //lammpsIO->lammpsdo(tolmp);
    
    //double *tempRad;
    //tempRad = new double[1];
    
    double tempRad = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"vidtt",0));
    
    //tempRad = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempRAD",1,1));
    //tempRad = ((double *) lammps_extract_atom(lammpsIO->lmp,(char *)"radius"));

    
    partVol = 4./3. * M_PI * tempRad*tempRad*tempRad;  // particle volume  FOR SPHERICAL

    
    //tolmp = "uncompute tempRAD";
    //lammpsIO->lammpsdo(tolmp);
    
    lammpsIO->lammpsdo("variable vidtt delete");
    lammpsIO->lammpsdo("uncompute craddd");
    lammpsIO->lammpsdo("uncompute tempRAD");
    
    
    tolmp = "delete_atoms group sub compress no";   // compress no avoids ids being re-assigned: important when fix_delete used in combination with fix_nucleate, so keep it!
    lammpsIO->lammpsdo(tolmp);
    tolmp = "group sub delete";
    lammpsIO->lammpsdo(tolmp);
    //tolmp = "run 1";   // just to see something in dump... later to be replaced with invoking a fix_relax
    //lammpsIO->lammpsdo(tolmp);
    
    lammpsIO->lammpsdo(lammpsIO->lmpThSt);    //restore original lammps thermo

    lammpsIO->lammpsdo( "undump tdID");
    
    
    
    
   
  
    

    //delete [] tempRad;
    
    
    
}



// ---------------------------------------------------------------
//  All procs (although only submaster is important...) adds elements to cumulative vector of number of particles in fix, and associated absolute ID of fix
void Fix_delete::add_Ng_apos(int Np, int apos)
{
    if (cumNPfix.size()==0)cumNPfix.push_back(Np);
    else cumNPfix.push_back(cumNPfix.back()+Np);
    
    fixaID.push_back(apos);
}



// ---------------------------------------------------------------
// Each submaster builds vectors with particle IDs, energies, rates etc, from all procs in subcomm. This function clears them
void Fix_delete::clearvecs()
{
    rates.clear();
    IDsrt.clear();
    ratesT.clear();
    CratesT.clear();
    cumNPfix.clear();
    fixaID.clear();
}


    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Fix_delete::printall()
{
	fprintf(screen,"\n---------ALL ABOUT FIX_DELETE----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
