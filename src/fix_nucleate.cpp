#include "fix_nucleate.h"
#include "fix.h"
#include "universe.h"
#include "chemistry.h"
#include "error.h"
#include "lammpsIO.h"
#include "output.h"
#include "krun.h"
#include "solution.h"
#include "store.h"
#include <memory.h>
#include <string.h>

using namespace MASKE_NS;

Fix_nucleate::Fix_nucleate(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    EVpID = -1;
    EVpXpos = -1.;
    EVpYpos = -1.;
    EVpZpos = -1.;
    EVpDIAM = -1.;
    EVpTYPE = -1;

    SAR = nullptr;
    Dsub = nullptr;
    CFuns = nullptr;
    GMuns = nullptr;
    fGMuns = nullptr;
    IDaruns = nullptr;
    Raruns = nullptr;

    //nmax = 0;
    //nlocal0 = 0;
    //x0 = NULL;
    
    id0 = new int[1];   // dummy memory allocation; will be deleted every time this array is populated anew
    typ0 = new int[1];   // dummy memory allocation; will be deleted every time this array is populated anew
    xyz0 = new double[1];   // dummy memory allocation; will be deleted every time this array is populated anew
}



// ---------------------------------------------------------------
// Class destructor
Fix_nucleate::~Fix_nucleate()
{
    //lammpsIO->lmp->memory->destroy(x0);
    delete [] id0;
    delete [] typ0;
    delete [] xyz0;
}







// --------------------------------------------------------------
void Fix_nucleate::delete_trial(int pos)
{
    std::string tolmp;
    
    // deleting trial particles in this fix
    std::ostringstream typ;
    typ << fix->fKMCptypeTRY[pos];
    tolmp = "group ptodel type "+typ.str();
    lammpsIO->lammpsdo(tolmp);
    tolmp = "delete_atoms group ptodel compress no";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "group ptodel delete";
    lammpsIO->lammpsdo(tolmp);
}



// ---------------------------------------------------------------
// Initialise the fix, see decription in fix_nucleate.h
void Fix_nucleate::init(int pos)
{
    // NB: lattice and region will not be evaluated until next KMC event is accepted. If meanwhile the box changes shape, the position of the trial particles should be updated automatically. This is best done making sure that trial particles are moved onto lattice just before calling a box-deforming event, and that any box deform applies an affine transformation to all particles. Also, this requires recording initial lattice cell volume and updating it after each box deform (if any). Best using the quickmin max distance parameters which can enable partly overlapping lattices to be devised. Then consider effect of deformation on volume dV/V = Tr(eps)
    std::string tolmp;
    std::ostringstream typ;
    
    // deleting trial particles in this fix
    /*
    tolmp = "group ptodel type "+typ.str();
    lammpsIO->lammpsdo(tolmp);
    tolmp = "delete_atoms group ptodel compress no";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "group ptodel delete";
    lammpsIO->lammpsdo(tolmp);
    */
    
    typ << fix->fKMCptypeTRY[pos];
    
    // activating fix-specific region and lattice in lammps
    int rid = fix->fKMCrid[pos];
    tolmp = store->RegCmd[rid];
    lammpsIO->lammpsdo(tolmp);
    int lid = fix->fKMClid[pos];
    tolmp = store->LatCmd[lid];
    lammpsIO->lammpsdo(tolmp);
    
    
    // creating particles of specified trial type on lattice in region (incl geometrical properties)
    tolmp = "create_atoms "+typ.str()+" region "+store->RegNames[rid];
    lammpsIO->lammpsdo(tolmp);
    if (strcmp((fix->fKMCpgeom[pos]).c_str(),"sphere")==0) {
        std::ostringstream diamt;
        diamt << fix->fKMCpdiam[pos];
        tolmp = "set type "+typ.str()+" diameter "+diamt.str();
        lammpsIO->lammpsdo(tolmp);
    } 
    
    
    // deactivating region in lammps (lattice cannot be deactivated, but is ok to be overwritten later)
    tolmp = "region "+store->RegNames[rid]+" delete";
    lammpsIO->lammpsdo(tolmp);

    
    // count the number of trial particles inserted = number of possible nucleate events in this fix
    // This will be useful later on, to decide where to count if an event in this fix is to be considered
    tolmp = "group gtemp type "+typ.str();
    lammpsIO->lammpsdo(tolmp);
    std::ostringstream stp;
    stp << msk->step;
    tolmp = "displace_atoms gtemp random 1E-8 1E-8 1E-8 "+stp.str()+" units box";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable tempNgroup equal count(gtemp)";
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");
    double tNg = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"tempNgroup",0));
    int Ng = (int) tNg;
    lammpsIO->lammpsdo("group gtemp delete");
    lammpsIO->lammpsdo("variable tempNgroup delete");
    fix->fKMCnevents[pos] = Ng;
    
    // computed DV associate to a particle on the lattice, via user-defined LAMMPS variable
    tolmp = "variable tempDV "+store->LatDVcmd[lid];
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");
    double tDV = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"tempDV",0));
    lammpsIO->lammpsdo("variable tempDV delete");
    fix->fKMC_DV[pos] = tDV;
    
    
    // exclude interaction of trial types in this fix with anything else (to be reactivated then just before sampling each fix nucleate, and deactivated immediately afterwards)
    for (int j=0; j<msk->Rtypes.size(); j++) {
        std::ostringstream typ2;
        typ2 << msk->Rtypes[j];
        tolmp = "neigh_modify exclude type "+typ.str()+" "+typ2.str();
        lammpsIO->lammpsdo(tolmp);
    }
    for (int j=0; j<msk->Ttypes.size(); j++) {
        std::ostringstream typ2;
        typ2 << msk->Ttypes[j];
        tolmp = "neigh_modify exclude type "+typ.str()+" "+typ2.str();
        lammpsIO->lammpsdo(tolmp);
    }
    typ.str("");
    typ.clear();

    
    // Record ID, type and x,y,z of all particles (target being only the trial, but faster to take all and filter out the trial later)
    if (pos == fix->fKMClast_nuc){  // if this is the last fix nucleate in this subcomm, then all trials have been created already
        
        
        // lammps extract all particles IDs, types and, via atom variables, x,y,z positions. Either use a collective lammps extraction, if you manage to get it to work, or do it per processor and then MPI_allgatherv to create a total array in all procs in this subcomm. Later, the reset trial pos function will go through these entire arrays and, if the ith atom is trial, it will use the id and x,y,z info to set position using LAMMPS interface.
        
        natoms0 = static_cast<int> (lammpsIO->lmp->atom->natoms); // total number of atoms in lammps (real and trial)
        
        delete [] id0;
        delete [] typ0;
        delete [] xyz0;
        id0 = new int[natoms0];
        typ0 = new int[natoms0];
        xyz0 = new double[natoms0*3];
        
        lammps_gather_atoms_concat(lammpsIO->lmp,(char *) "id",0,1,id0);
        lammps_gather_atoms_concat(lammpsIO->lmp,(char *) "type",0,1,typ0);
        lammps_gather_atoms_concat(lammpsIO->lmp,(char *) "x",1,3,xyz0);

        if (msk->wplog) {
            std::string msg = "\nAll particles data from latest fix_nucleate.init() in this subcomm";
            output->toplog(msg);
            std::ostringstream ss;
            ss <<"nTotal number of particles (real + all trial): "<<natoms0;
            msg = ss.str();     ss.str("");     ss.clear();
            output->toplog(msg);
            msg = "pos id type x y z";
            output->toplog(msg);
            for (int i=0; i<natoms0; i++){
                ss << i<<" "<< id0[i] <<" "<<typ0[i]<<" "<<xyz0[3*i]<<" "<<xyz0[3*i+1]<<" "<<xyz0[3*i+2];
                msg = ss.str();     ss.str("");     ss.clear();
                output->toplog(msg);
            }
        }
    }
    
    
    // record vectors with initial x,y,z positions of all particles but not yet the initial box size. This info will be used at each call of sample later on to apply affine deformations and adjust volume per trial particle until init is invoked again
    /*
    if (nmax < lammpsIO->lmp->atom->nmax) {
      tag0 = lammpsIO->lmp->memory->grow(tag0,lammpsIO->lmp->atom->nmax,"fix_nucleate:tag0");
      type0 = lammpsIO->lmp->memory->grow(type0,lammpsIO->lmp->atom->nmax,"fix_nucleate:type0");
      x0 = lammpsIO->lmp->memory->grow(x0,lammpsIO->lmp->atom->nmax,3,"fix_nucleate:x0");
    }
    nlocal0 = lammpsIO->lmp->atom->nlocal;
    // assuming continuous chunck of memory allocated through memory->grow 
    memcpy(tag0,lammpsIO->lmp->atom->tag,nlocal0*sizeof(LAMMPS_NS::tagint));
    memcpy(type0,lammpsIO->lmp->atom->type,nlocal0*sizeof(int));
    memcpy(x0[0],lammpsIO->lmp->atom->x[0],nlocal0*3*sizeof(double));
     */
}

// ---------------------------------------------------------------
// Reset trial atoms initial positions. This assumes that there was no nucleation or dissolution happening between call to init() and reset_trial_pos(), hence we can still use atom->nlocal as the number of atoms for the stored positions.
void Fix_nucleate::reset_trial_pos(int pos)
{
    std::ostringstream ss;
    for (int i=0; i<natoms0; i++){
        if (typ0[i] == fix->fKMCptypeTRY[pos]) {
            ss <<"set atom "<<id0[i]<<" x "<<xyz0[3*i]<<" y "<<xyz0[3*i+1]<<" z "<<xyz0[3*i+2];
            lammpsIO->lammpsdo(ss.str());
            ss.str("");     ss.clear();
        }
    }
    
    /*
    int nwarn = 0;
    for (int i = 0; i < nlocal0; i++) {
        if (type0[i] == fix->fKMCptypeTRY[pos]) {
          int index = lammpsIO->lmp->atom->map(tag0[i]);
          if (index >= 0) {
            lammpsIO->lmp->atom->x[index][0] = x0[i][0];
            lammpsIO->lmp->atom->x[index][1] = x0[i][1];
            lammpsIO->lmp->atom->x[index][2] = x0[i][2];
          } else if(nwarn==0){
              fprintf(screen, "[WARNING] Can't find atom %d in Fix_nucleate::reset_trial_pos on rank %d",tag0[i],universe->me);
              nwarn++;
          }
        }
    }
     */
}

// ---------------------------------------------------------------
// Sampling all delete transitions correspondin to fix in position "pos" within this subcommunicator list, recorded in the fix object
void Fix_nucleate::sample(int pos)
{
    // USEFUL SHORTCUTS for current processor key and number of procs in local subcomm
    key = universe->key;
    nploc = universe->SCnp[universe->color];
    
    mid = fix->fKMCmid[pos];    // mechanism ID
    
    std::string tolmp;
    
    // screen output to debug
    int sleeptime = me;
    if (msk->wplog) {
        std::string msg = "\nTEMPO ";
        std::ostringstream ss;
        ss << msk->tempo;
        msg = msg+ss.str()+": sampling fix_nucleate\""+fix->fKMCname[pos]+"\"";
        output->toplog(msg);
    }
    
    
    
    
    // TRIAL PARTICLES SEARCH FOR LOCAL MINIMUM OF INTERACTION ENERGY
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    
    // create temp group of all real particles, and freeze them
    tolmp = "group gtempFIX type ";
    for (int i=0; i<msk->Rtypes.size(); i++) {
        std::ostringstream typ;
        typ << msk->Rtypes[i];
        tolmp = tolmp + typ.str() + " ";
    }
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "fix tmpfreeze gtempFIX setforce 0.0 0.0 0.0";
    //tolmp = "fix tmpfreeze gtemp freeze";
    lammpsIO->lammpsdo(tolmp);
    
    
    // exclude interactions
    tolmp = "neigh_modify exclude none";   // start with all particles seeing all, then subtract..
    lammpsIO->lammpsdo(tolmp);
    
    // exclude interaction of trial particles in this fix among themselves
    std::ostringstream typ;
    typ << fix->fKMCptypeTRY[pos];
    tolmp = "neigh_modify exclude type "+typ.str()+" "+typ.str();
    lammpsIO->lammpsdo(tolmp);
    typ.str("");
    typ.clear();
    
    // exclude interactions of all other trial particles except for those in this fix
    for (int i=0; i<msk->Ttypes.size(); i++) {
        if (msk->Ttypes[i] != fix->fKMCptypeTRY[pos]) {
            typ << msk->Ttypes[i];
            for (int j=0; j<msk->Rtypes.size(); j++) {
                std::ostringstream typ2;
                typ2 << msk->Rtypes[j];
                tolmp = "neigh_modify exclude type "+typ.str()+" "+typ2.str();
                lammpsIO->lammpsdo(tolmp);
            }
            for (int j=0; j<msk->Ttypes.size(); j++) {
                std::ostringstream typ2;
                typ2 << msk->Ttypes[j];
                tolmp = "neigh_modify exclude type "+typ.str()+" "+typ2.str();
                lammpsIO->lammpsdo(tolmp);
            }
            typ.str("");
            typ.clear();
        }
    }
    
    
    // minimise
    int minid = fix->fKMCminid[pos];
    tolmp = "timestep "+store->MinTstep[minid];
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "min_style quickmaske";
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = store->MinModCmd[minid];
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = store->MinCmd[minid];
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "timestep 0";
    lammpsIO->lammpsdo(tolmp);
    
    // END OF POSITIONING TRIAL PARTICLES IN LOCAL MINIMA OF ENERGY
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    
    
    
    
    


    // EXTRACTING IDs AND ENERGIES FROM LAMMPS
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    
    
    tolmp = "compute tempID all property/atom id"; // temp compute to get id of all atoms
    lammpsIO->lammpsdo(tolmp);
    
    // radii are needed for nucleation onle for coverage fraction in micro mechanism: you can probable improve efficiency by doing this only if a micro mechanism is invoked and not in general (NB: in fix delete instead you need radii to compute number of layers in general, even for non micro mechanisms)
    tolmp = "compute tempRAD all property/atom radius"; // temp compute to get atom radii
    lammpsIO->lammpsdo(tolmp);
    

    typ << fix->fKMCptypeTRY[pos];
    tolmp = "group gtemp type "+typ.str();
    typ.str(""); typ.clear();
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable tempNgroup equal count(gtemp)";
    //tolmp = "variable tempNgroup equal count("+fix->fKMCgroup[pos]+")";  // temp variable counting atoms of type in fix
    lammpsIO->lammpsdo(tolmp);
    
    // NB: IF POTENTIAL IS NOT PAIRWISE, but e.g. three body, will need to change relationship between pe/atom extracted here and DU computed later for the rate. If potential is a mix of pair and other terms, e.g. a generic one for molecular interactions, I will need to compute energies per particle or molecule in a different way, by first minimising all trial particles, then assiging all but one trial to a group with deactivated interactions, computing the DUtot with only one active trial, and then looping over all trials - more expensive
    tolmp = "compute tempPE all pe/atom";  // temp compute reading energy per atom
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "compute tempType all property/atom type";  // temp compute reading type per atom
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "dump tdID all custom 1 dump.tempTRY_"+universe->SCnames[universe->color]+" c_tempID c_tempPE type x y z";    //temp dump to update variables and computes -- NOTICE I removed c_tempRAD compared to fix_delete
    lammpsIO->lammpsdo(tolmp);
    
    if(strcmp((chem->mechstyle[mid]).c_str(),"micro")==0){
        if(strcmp((chem->mechpar[mid][0]).c_str(),"pair")==0){
            tolmp = "compute tempPAT all property/local patom1 patom2 ptype1 ptype2";
            lammpsIO->lammpsdo(tolmp);
            tolmp = "compute tempDIST all pair/local dist";
            lammpsIO->lammpsdo(tolmp);
            tolmp = " dump tdLID all local 1 dumpL.tmp index c_tempPAT[1] c_tempPAT[2] c_tempPAT[3] c_tempPAT[4] c_tempDIST";
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
        std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: fix_nucleate.cpp found "+wd4.str()+" atoms in fix "+fix->fKMCname[pos]+", whereas krun.cpp previously found "+wd5.str()+" atoms in same fix. This should never happen. Something is wrong with the source code... \n****************************************************\n***********************************************\n";
        fprintf(screen,"%s",msg.c_str());
        //error->errsimple(msg);
    }
    
   
     
    natoms = static_cast<int> (lammpsIO->lmp->atom->natoms); // total number of atoms in lammps, all groups all fixes)
    nlocal = static_cast<int> (lammpsIO->lmp->atom->nlocal); // number of atoms in current processor, all types all fixes)
    
    //aID = new double[natoms];     // vector of unsorted IDs
    aID = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempID",1,1));
    // NB: each processor in subcom pulls out the ID of their atoms. We will put them all into a single vector, IDuns, to be managed by the submaster. The way that seems to work is to scan aID of each processor looking for the first nlocal atoms with non-zero id. Ids after those are random. The first nonzero nlocal ids are passed to the submaster, which eventually sorts them.
    
    // arrays of radii and energies
    //aR = new double[natoms];
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
            #ifdef MASKE_WITH_NUFEB 
            nlocR = *((int*)lammps_extract_compute(lammpsIO->lmp,(char *) "tempPAT",2,0));
            #else
            nlocR = *((int*)lammps_extract_compute(lammpsIO->lmp,(char *) "tempPAT",2,4));
            #endif
        }
    }

    
    
    //printing group-specific aIDs and all IDs in lammps vectors, for debugging
    if (msk->wplog) {
        std::string msg = "Total number of atoms in LAMMPS: ";
        std::ostringstream ss;
        ss << natoms;       msg = msg+ss.str(); ss.str("");   ss.clear();
        msg = msg+"\nNumber of atoms in current processor: ";
        ss << nlocal;       msg = msg+ss.str(); ss.str("");   ss.clear();
        output->toplog(msg);
        msg="";
        output->toplog("\npos atype aiD aR aE");
        {
            int naf = 0;  // counter of number of atoms in current processor
            int i = 0;
            while (naf < nlocal){
                if (aID[i]>0){
                    naf++;
                    if ((int)atype[i]==fix->fKMCptypeTRY[pos]){
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
    

    tolmp = "undump tdID ";     // removing temporary dump
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
    
    
    
    
    

    
    
    // RECORDING IDs AND ENERGIES IN EACH PROCESSOR (cleaning out those not pertaining to current processor)
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    tID.clear();    tE.clear();   tR.clear();
    naP = 0;  // number of atoms in current processor of type invoked by current fix
    int naf = 0;  // counter of number of atoms in current processor
    {
        int i = 0;
        while (naf < nlocal){
            if (aID[i]>0){
                naf++;
                if ((int)atype[i]==fix->fKMCptypeTRY[pos]){
                    tID.push_back((int)aID[i]);
                    tR.push_back(aR[i]);
                    tE.push_back(aE[i]);
                    naP++;
                }
            }
            i++;
        }
    }

    
    
    //delete [] aID;

    //delete [] aE;
    //delete [] atype;

    
    tIDarr = new int[naP];     // copying local ID vector to array in order to communicate it to submaster later on
    for (int i =0; i<naP; i++) tIDarr[i] = tID[i];
    
    //printing proc-specific tIDs (both vecto and array), and all IDs in lammps vectors, for debugging
    if (msk->wplog) {
        std::string msg = "\nNumber of atoms in this processor: ";
        std::ostringstream ss;    ss << naP;   msg = msg+ss.str();
        output->toplog(msg);
        msg="";   ss.str("");   ss.clear();
        output->toplog("\npos tiD tiDarr tE");
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
                std::string msg = "\nNumber of real atoms in this processor (for coverage function in fix_nucleate.cpp): ";
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
    else if (strcmp((chem->mechstyle[mid]).c_str(),"micro")==0) {
        comp_rates_micro(pos);
    }
    //else if (strcmp((chem->mechstyle[mid]).c_str(),"allser")==0) {
      //  comp_rates_allser(pos);
    //}
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
    
    // release frozen real particles
    tolmp = "unfix tmpfreeze";
    lammpsIO->lammpsdo(tolmp);
    
    // exclude interaction of trial types in this fix with anything else
    typ << fix->fKMCptypeTRY[pos];
    for (int j=0; j<msk->Rtypes.size(); j++) {
        std::ostringstream typ2;
        typ2 << msk->Rtypes[j];
        tolmp = "neigh_modify exclude type "+typ.str()+" "+typ2.str();
        lammpsIO->lammpsdo(tolmp);
    }
    for (int j=0; j<msk->Ttypes.size(); j++) {
        std::ostringstream typ2;
        typ2 << msk->Ttypes[j];
        tolmp = "neigh_modify exclude type "+typ.str()+" "+typ2.str();
        lammpsIO->lammpsdo(tolmp);
    }
    typ.str("");
    typ.clear();

    lammpsIO->lammpsdo("group gtempFIX delete");
    
    
    delete [] rate_each;
    delete [] tIDarr;
    delete [] nID_each;
    delete [] IDpos;
    delete [] IDuns;
    
    
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

            if (!(Dsub==nullptr)){
                delete [] Dsub;
                Dsub=nullptr;
            }
            if (!(CFuns==nullptr)){
                delete [] CFuns;
                CFuns=nullptr;
            }
            if (!(GMuns==nullptr)){
                delete [] GMuns;
                GMuns=nullptr;
            }
            if (!(fGMuns==nullptr)){
                delete [] fGMuns;
                fGMuns=nullptr;
            }
            if (!(IDaruns==nullptr)){
                delete [] IDaruns;
                IDaruns=nullptr;
            }
            if (!(Raruns==nullptr)){
                delete [] Raruns;
                Raruns=nullptr;
            }
            
            delete [] IDar;
            delete [] Rar;
            delete [] nIDar_each;
            delete [] IDarpos;
            
        }
    }

    
    
}




// ---------------------------------------------------------------
// assembles all IDs from each processor into block-wise array in submaster, where each block is the tIDarr in a processor
void Fix_nucleate::ids_to_submaster(int pos)
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
    //fprintf(screen,"\n\n  PROC %d , SUBCOM %d , fix %s, tempo %f\n I have %d local ids   \n",me,universe->color,fix->fKMCname[pos].c_str(),msk->tempo,nID_each[key]);
    
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
    IDarpos = new int[nploc];  // position of local real particle array in submaster's unsorted list of all real IDs

    if (key==0) {
        //fprintf(screen,"\n\n  PROC %d , SUBCOM %d , fix %s, tempo %f\n Here is the nID_each that I know of:  \n",me,universe->color,fix->fKMCname[pos].c_str(),msk->tempo);
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
void Fix_nucleate::pair_arr_to_submaster(int pos)
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
void Fix_nucleate::submaster_sort_IDs(int pos)
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
void Fix_nucleate::submaster_map_ID(int pos)
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
        output->toplog("\nIDsrt");
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
void Fix_nucleate::submaster_comp_cover(int pos)
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
        output->toplog("npos id1 id2 type1 type2 upos1 upos2 Dsub");
    }
    
    //printf(screen,"DEBUG 1: PROC %d try type is %d",me,fix->fKMCptypeTRY[pos]);
    //sleep (1);
    
    
    // for all the rows in the assembled arrays in the submaster...
    for (int i=0; i<SARpos[nploc-1]+nlocR_each[nploc-1]; i++){
        
        //sleep(1);
        //std::cout << "\n DEBUG 1 " << "proc " << me << std::endl;
        //sleep(1);
        
        int id1 = SAR[i][0];
        int id2 = SAR[i][1];
        int t1 = SAR[i][2];
        int t2 = SAR[i][3];
        bool flag_t1 = (t1 == fix->fKMCptypeTRY[pos]);
        bool flag_t2 = (t2 == fix->fKMCptypeTRY[pos]);

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
            ss << Dsub[i];   msg += ss.str()+" ";  ss.str("");   ss.clear();
            output->toplog(msg);
        }
        
        
      
        if (flag_t1 || flag_t2){
            // compute coverage area, used later to compute number of layers to dissolve in micro particle
            // coverage is given by contact cross section weighted by a distance-dependent factor; the latter is user-provided through two thresholds, e0 below which contact is 100%, ef above which contact is 0%. Linear interpolation in between
            
            double Aij;    // cross section of inter-particle contact
            double Rij;     // harmonic average of radii of particles in contact
            double Ri,Rj;   // radii of interaction particles in current pair
            int rp1,rp2;    // position of interacting particles in the "all real" sorted vectors (if not of fix-sampled trial type, for which radius is assigned by fix itself)
            
            //sleep(1);
            //std::cout << "\n DEBUG 2 " << "proc " << me << std::endl;
            //sleep(1);
            
            // find first particle position in sorted IDar vector
            if (flag_t1) Ri = fix->fKMCpdiam[pos]/2.;
            else{
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
                Ri = Rarsrt[rp1];
            }
          
            //sleep(1);
           // std::cout << "\n DEBUG 3 " << "proc " << me << std::endl;
            //sleep(1);
            
            if (flag_t2) Rj = fix->fKMCpdiam[pos]/2.;
            else{
                // find second particle position in sorted IDar vector
                int pre = 0;
                int post = IDarsrt.size()-1;
                int posit;
                while (pre < post) {
                    posit = (int)((pre+post)/2.);
                    if (IDarsrt[posit] < id2) pre=posit+1;
                        else post=posit;
                    }
                posit=pre;
                rp2 = posit;
                Rj = Rarsrt[rp2];
            }
            
            
            Rij = 2. * Ri * Rj/( Ri + Rj);
            Aij = M_PI * Rij * Rij;
            
            //sleep(1);
            //std::cout << "\n DEBUG 4 " << "proc " << me << std::endl;
            //sleep(1);
            
            // weigh the contact cross section by the distance
            double Dij; // arithmetic average of diameters in contact
            Dij = Ri + Rj;
            double efij = chem->ef[t1-1][t2-1];
            double e0ij = chem->e0[t1-1][t2-1];
            
            if (Dsub[i] > efij * Dij)        Aij = 0.;
            else if (Dsub[i] > e0ij * Dij)   Aij *= (efij - Dsub[i]/Dij) / (efij- e0ij);

           // sleep(1);
           // std::cout << "\n DEBUG 5 " << "proc " << me << std::endl;
           // sleep(1);
           
            // if first atom is of correct type for this fix, add contact area to the contact fraction arrays (still areas here; will be converted to area fractions later on below)
            if (flag_t1)  {
                CFuns[up1] += Aij / (4. * M_PI * Ri * Ri);
                if (Aij > 0.){
                    double tempGM = 0.;
                    if (t1 != t2) tempGM =  chem->gij[t1-1][t1-1] + chem->gij[t1-1][t2-1] - chem->gij[t2-1][t2-1];
                    if (fGMuns[up1]) {
                        GMuns[up1] = tempGM;
                        fGMuns[up1] = false;
                    }
                    else if (tempGM > GMuns[up1]) GMuns[up1] = tempGM;
                }
            }
            
           // sleep(1);
            //std::cout << "\n DEBUG 6 " << "proc " << me << std::endl;
            //sleep(1);
            
            if (flag_t2){
                CFuns[up2] += Aij / (4. * M_PI * Rj * Rj);
                if (Aij > 0.){
                    double tempGM = 0.;
                    if (t1 != t2) tempGM =  chem->gij[t2-1][t2-1] + chem->gij[t2-1][t1-1] - chem->gij[t1-1][t1-1];
                    if (fGMuns[up2]) {
                        GMuns[up2] = tempGM;
                        fGMuns[up2] = false;
                    }
                    else if (tempGM > GMuns[up2]) GMuns[up2] = tempGM;
                }
            }
            
            //std::cout << "\n DEBUG 7 " << "proc " << me << std::endl;
        }
        
    }
    double nbulk = 12.;   // number of particles in the bulk for a monodisperse system of particls with type selected for dissolution. In future implementations, this could be passed as an input when calling the fix_delete (becasue it applies to all particles to be considered for deletion here)
    for (int i=0; i<Ng; i++){
        // convert coverage areas to coverage fractions
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
void Fix_nucleate::cover_from_submaster(int pos)
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
void Fix_nucleate::comp_rates_allpar(int pos)
{
    
    // for all particles in current ("each") processor
    for (int i =0; i<nID_each[key]; i++) {
        
        if (msk->wplog) {
            std::string msg = "\nPARTICLE ";
            std::ostringstream ss;    ss << tID[i];   msg = msg+ss.str();
            output->toplog(msg);
        }
        
        double Pv, Ps;  // particle volume and surface
        if (strcmp(fix->fKMCpgeom[pos].c_str(),"sphere")==0) {
            double Pd = fix->fKMCpdiam[pos];
            Pv = M_PI / 6. * Pd * Pd * Pd;
            Ps = M_PI * Pd * Pd;
        }

        int nrt = 1;    //number of reaction in series in chain
        double nrv = 1;    // number of unit reactions or unit chains in particle volume
        
        int chID = -1;  // chain ID: will be > -1 only if mechanism calls indeed a chain
        int rxid = -1;  // reaction ID: will be used now but later in loop too
        
        bool errflag = false;   // error to spot if wrong volume changes defined for reactions: see below
        double fgdV = 0;
        if (chem->mechchain[mid]) {  //if reaction is a chain
            chID = chem -> mechrcID[mid];
            nrt = (chem->ch_rxID[chID]).size();
            if (chem -> ch_dVp_fgd[chID] > 0) {
                nrv = round( Pv / chem -> ch_dVp_fgd[chID]);   // NB: a proper nucleation reax should have positive fgd changes defined in the chemDB file
            } else {errflag = true; fgdV=chem -> ch_dVp_fgd[chID];}
        }
        else{ //if instead it is a single reaction
            nrt = 1;
            rxid = chem->mechrcID[mid];
            if (chem -> rx_dVp_fgd[rxid] > 0) {
                nrv = round( Pv / chem -> rx_dVp_fgd[rxid]);
            } else {errflag = true; fgdV=chem -> rx_dVp_fgd[rxid];}
        }
        
        if (errflag == true) {
            std::ostringstream wd1,wd2,wd3,wd4;
            wd1<<me; wd2<<universe->color; wd3<<msk->tempo,wd4<<fgdV;
            std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: fix "+fix->fKMCname[pos]+" is of nucleation type, which requires positive fgd volume changes, whereas its associated mechanism has fgd volume change = "+wd4.str()+" \n****************************************************\n***********************************************\n";
            fprintf(screen,"%s",msg.c_str());
        }
        
        std::string msg;
        if (msk->wplog) { msg = "\n nrv "; std::ostringstream ss;    ss << nrv;   msg = msg+ss.str(); }
        
        // change of surface energy and interaction energy per reaction unit (single reaction or chain: to be further subdivided per step in chain later on)
        double DSpu;  // positive becasue nucleation
        double DUpu;;  //
        
        if (strcmp(chem->mechinter[mid].c_str(),"int_1lin")==0)  {
            DSpu = Ps/pow((double)nrv,1./3.);
            DUpu = 2.*tE[i]/pow((double)nrv,1./3.);
        }
        else if (strcmp(chem->mechinter[mid].c_str(),"int_2lin")==0) {
            DSpu = Ps/pow((double)nrv,68.5/100.); // power artificailly tuned to get hetero nucleation in cg sim.
            DUpu = 2.*tE[i]/pow((double)nrv,68.5/100.);
        }
        else{
            DSpu = Ps/((double)nrv);
            DUpu = 2.*tE[i]/((double)nrv);
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
                double gammax = chem -> compgammax(rxid);   // somewhere in this function I will need to calculate gammax from the concentration of the activated complex in case the complex is not charge neutral
                double KT = msk->kB * solution->Temp;
                                
                if (msk->wplog) {
                    std::ostringstream ss;
                    msg += "; DG* ";    ss << DGx;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; gamma* ";    ss << gammax;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; KB ";    ss << (msk->kB);   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; T ";    ss << (solution->Temp);   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                //double r0star =  KT / msk->hpl / gammax * exp(-DGx / KT);
                //double apf = 1.; // molecule size in prefactor: for "allpar" one should take 1 and use gammax above in dimensionless units.   PROBABLY USELESS unless using the Shvab's mechanism of nucleation by growth
                
                double kappa = 1.;  // transmission coefficient
                double Sen = chem->compSen(mid);
                
                
                //bool net = false; // the beta calculator in simulation.cpp uses stoichio coefficients that are defined by user in chemDB for forward reaction. If beta is computed for the net rate, the reaction is backward hence stoichio coeff need multiplication by -1   ... TO BE REVISED after new expressions of rates
                std::string topass = "reac";
                double Qreac = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]); // using activity products instead of supersaturation beta, which gives more generality. See Notes_on_TST.pdf document
                
                double Vt = (chem -> rx_dVt_fgd[rxid]);  // Volume of fgd created by the reaction at the current step in the chain. Plus sign because dV of fdg > 0 for a proper nucleation reaction.

                
                double DV =  fix->fKMC_DV[pos];
                
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
                    msg += "; DV ";    ss << DV;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Qreac ";    ss << Qreac;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Sen ";    ss << Sen;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; DSi ";    ss << DSi;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; DUi ";    ss << DUi;   msg += ss.str();    ss.str(""); ss.clear();
                    msg += "; Vti ";    ss << Vt;   msg += ss.str();    ss.str(""); ss.clear();
                }
                
                //double fa = 1.;
                //if (dim == 2 || dim== 1) fa = 3;
                
                r0 = r0*pow(Vt,(dim/3. - 1.))*DV;
                
                double ki = (chem -> ki[rxid]);
                
                // Rate equation as per TST (see notes_on_TST.pdf by Masoero, 2019)
                ri = r0 * Qreac     *exp(ki*(-Sen * DSi - DUi)/KT) ;

                if (msk->wplog) {
                    msg += ", ri ";
                    std::ostringstream ss;    ss << ri;   msg += ss.str(); ss.str("");   ss.clear();
                }
                
                
                // if net rate is requested by user in chemDB, then add rate backward
                if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                    
                    
                    topass = "prod";
                    double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    
                    if (msk->wplog) {
                        msg += ", Qprod ";
                        std::ostringstream ss;    ss << Qprod;   msg += ss.str(); ss.str("");   ss.clear();
                    }
                    
                    //net = true;
                    //beta = solution->compbeta(rxid,net,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    ri -= r0 * Qprod / chem->Keq[rxid] * exp((1.-ki)*(Sen * DSi + DUi)/KT ) ;
                    
                }
                
                if (msk->wplog) {
                    msg += ", rinet ";
                    std::ostringstream ss;    ss << ri;   msg += ss.str(); ss.str("");   ss.clear();
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
            msg += " has nucleation rate  ";
            ss << rate_each[i];   msg += ss.str();
            output->toplog(msg);
        }
    }
    
}



// ---------------------------------------------------------------
//  each processor computes its deletion rates
void Fix_nucleate::comp_rates_micro(int pos)
{
    
    // for all particles in current ("each") processor
    for (int i =0; i<nID_each[key]; i++) {
        
        if (msk->wplog) {
            std::string msg = "\n -----------------------------------\nPARTICLE ";
            std::ostringstream ss;    ss << tID[i];     msg = msg+ss.str()+" ";
            ss.str("");     ss.clear();   ss << tR[i];  msg = msg+ss.str()+" ";
            ss.str("");     ss.clear();   ss << tCF[i];  msg = msg+ss.str()+" ";
            ss.str("");     ss.clear();   ss << tGM[i];  msg = msg+ss.str()+"\n";
            output->toplog(msg);
        }
        
        double Pv;  // particle volume
        if (strcmp(fix->fKMCpgeom[pos].c_str(),"sphere")==0) {
            double Pd = fix->fKMCpdiam[pos];
            Pv = M_PI / 6. * Pd * Pd * Pd;
        }

        int nrt = 1;    //number of reaction in series in chain
        double nrv = 1;    // number of unit reactions or unit chains in particle volume
        
        int chID = -1;  // chain ID: will be > -1 only if mechanism calls indeed a chain
        int rxid = -1;  // reaction ID: will be used now but later in loop too
        
        bool errflag = false;   // error to spot if wrong volume changes defined for reactions: see below
        double fgdV = 0;
        if (chem->mechchain[mid]) {  //if reaction is a chain
            chID = chem -> mechrcID[mid];
            nrt = (chem->ch_rxID[chID]).size();
            if (chem -> ch_dVp_fgd[chID] > 0) {
                nrv = round( Pv / chem -> ch_dVp_fgd[chID]);   // NB: a proper nucleation reax should have positive fgd changes defined in the chemDB file
            } else {errflag = true; fgdV=chem -> ch_dVp_fgd[chID];}
        }
        else{ //if instead it is a single reaction
            nrt = 1;
            rxid = chem->mechrcID[mid];
            if (chem -> rx_dVp_fgd[rxid] > 0) {
                nrv = round( Pv / chem -> rx_dVp_fgd[rxid]);
            } else {errflag = true; fgdV=chem -> rx_dVp_fgd[rxid];}
        }
        
        if (errflag == true) {
            std::ostringstream wd1,wd2,wd3,wd4;
            wd1<<me; wd2<<universe->color; wd3<<msk->tempo,wd4<<fgdV;
            std::string msg = "****************************************************\n***********************************************\n PROC "+wd1.str()+" SUBCOM "+wd2.str()+" tempo "+wd3.str()+" \n ERROR: fix "+fix->fKMCname[pos]+" is of nucleation type, which requires positive fgd volume changes, whereas its associated mechanism has fgd volume change = "+wd4.str()+" \n****************************************************\n***********************************************\n";
            fprintf(screen,"%s",msg.c_str());
        }
        
        
        // change of interaction energy per reaction unit (single reaction or chain: to be further subdivided per step in chain later on)
        double DUpu;;  //
        
        if (strcmp(chem->mechinter[mid].c_str(),"int_1lin")==0)  {
            //DSpu = Ps/pow((double)nrv,1./3.);
            //DUpu = 2.*tE[i]/pow((double)nrv,1./3.);
        }
        else if (strcmp(chem->mechinter[mid].c_str(),"int_2lin")==0) {
            //DSpu = Ps/pow((double)nrv,68.5/100.); // power artificailly tuned to get hetero nucleation in cg sim.
            //DUpu = 2.*tE[i]/pow((double)nrv,68.5/100.);
        }
        else{
            //DSpu = Ps/((double)nrv);
            DUpu = 2.*tE[i]/((double)nrv);
        }
        
        
        // number of layers to dissolve (it will depend on fractional coverage area)
        double lim_CF= std::stod(chem->mechpar[mid][3]); //limiting coverage fraction
        
        int nrL = 0;
        
        double unit_thick;    // thickness of precipitating unit in radial direction of particle
        if (chem->mechchain[mid]) {  //if reaction is a chain
            unit_thick = pow(chem->ch_dVp_fgd[chID],1./3.);
        }
        else{ //if instead it is a single reaction
            unit_thick = pow(chem->rx_dVp_fgd[rxid],1./3.);
        }
        
        if (tCF[i] > 0.5) nrL = round(tR[i]/unit_thick);
        else if (tCF[i] >= lim_CF) nrL = round(2.*tR[i]/unit_thick);
        
        bool flag_bulk = false;
        if (nrL == 0) flag_bulk = true;     // if particle has < lim_CF coverage fraction it is nucleating homogeneously in the bulk solution: not allowed in this mechanism. This scenario should be sampled separately with another dedicated mechanisms, e.g. CNT
        
        
        // number of units to precipitate a full layer (depends on density of kinks per layer)
        int nrS = 1;
        if (chem->mechchain[mid])  nrS = round(1./chem->ch_Fk[chID]);
        else nrS = round(1./chem->rx_Fk[rxid]);
        
        
        
        // compute rates for all reaction units (multi-step if a chain), for all layers, for all surface stpes
        
        double ri=0. , DTi = 0., DTtot = 0.;  // rate of each reaction in sequence, associated time increement and total cumulative time increment
        
        //=============================================MICRO LAYER-BY-LAYER SLOW RUN (START)===========================================================
        /*for (int j=0; j<nrL-1; j++){
            
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
                    
                    if (msk->wplog) { msg += "; rx_step ";    std::ostringstream ss;    ss << k;   msg += ss.str(); }
                    
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
                    
                   
                    double kappa = 1.;  // transmission coefficient
                    
                    std::string topass = "reac";
                    double Qreac = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    
                    double Vt = (chem -> rx_dVt_fgd[rxid]);  // Tributary volume of fgd deleted by the reaction at the current step in the chain.

                    // lattice cell volume: needed for rate
                    double DV =  fix->fKMC_DV[pos];
                    
                    
                    double DUi = 0.;
                    if (strcmp(chem->mechinter[mid].c_str(),"int_no")!=0) {
                        DUi = DUpu;
                        if (chem->mechchain[mid]) DUi *= (chem->ch_rdV_fgd[chID][k]);
                        
                        // Add change in strain energy internal to the molecules
                        for (int j1=0; j1<chem->fgd_molID[rxid].size(); j1++){
                            int molid = chem->fgd_molID[rxid][j1];
                            DUi += chem->mol_Us[molid] * chem->mol_vm[molid] * chem->fgd_nmol[rxid][j1];
                        }
                    }
                    
                    
                    
                    
                    double r0 = kappa * KT / msk->hpl / gammax * cx * exp(-DGx / KT);
                    
                    if (msk->wplog) {
                        std::ostringstream ss;
                        msg += "; r0 ";    ss << r0;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; DV ";    ss << DV;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Qreac ";    ss << Qreac;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; DUi ";    ss << DUi;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Vti ";    ss << Vt;   msg += ss.str();    ss.str(""); ss.clear();
                    }
                    
                    r0 = r0*pow(Vt,(dim/3. - 1.))*DV;
                    
                    double ki = (chem -> ki[rxid]);
                    
                    // Forward rate equation
                    ri = r0 * Qreac * exp(ki * (- DUi)/KT ) ;
                    

                    
                    // if net rate is requested by user in chemDB, then add rate backward
                    if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                        
                        
                        topass = "prod";
                        double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                        
                        ri -= r0 * Qprod / chem->Keq[rxid] * exp((1.-ki)*(DUi)/KT ) ;
                        
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
            
            if (msk->wplog) {
                msg += "\n";
                output->toplog(msg);
            }
        }*/
        //=============================================MICRO LAYER-BY-LAYER SLOW RUN (END)===========================================================

        //=====================================================MICRO FAST RUN (START)==================================================================

        // compute rate of reaction sequence
                for (int k=0; k< nrt; k++) { //all the reaction in series in chain seq.
                    std::string msg;
                    
                    if (msk->wplog) { msg += "; rx_step ";    std::ostringstream ss;    ss << k;   msg += ss.str(); }
                    
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
                    
                   
                    double kappa = 1.;  // transmission coefficient
                    
                    std::string topass = "reac";
                    double Qreac = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    
                    double Vt = (chem -> rx_dVt_fgd[rxid]);  // Tributary volume of fgd deleted by the reaction at the current step in the chain.

                    // lattice cell volume: needed for rate
                    double DV =  fix->fKMC_DV[pos];
                    
                    
                    double DUi = 0.;
                    if (strcmp(chem->mechinter[mid].c_str(),"int_no")!=0) {
                        DUi = DUpu;
                        if (chem->mechchain[mid]) DUi *= (chem->ch_rdV_fgd[chID][k]);
                        
                        // Add change in strain energy internal to the molecules
                        for (int j1=0; j1<chem->fgd_molID[rxid].size(); j1++){
                            int molid = chem->fgd_molID[rxid][j1];
                            DUi += chem->mol_Us[molid] * chem->mol_vm[molid] * chem->fgd_nmol[rxid][j1];
                        }
                    }
                    
                    
                    
                    
                    double r0 = kappa * KT / msk->hpl / gammax * cx * exp(-DGx / KT);
                    
                    if (msk->wplog) {
                        std::ostringstream ss;
                        msg += "; r0 ";    ss << r0;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; DV ";    ss << DV;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Qreac ";    ss << Qreac;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; DUi ";    ss << DUi;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Vti ";    ss << Vt;   msg += ss.str();    ss.str(""); ss.clear();
                    }
                    
                    r0 = r0*pow(Vt,(dim/3. - 1.))*DV;
                    
                    double ki = (chem -> ki[rxid]);
                    
                    // Forward rate equation
                    ri = r0 * Qreac * exp(ki * (- DUi)/KT ) ;
                    

                    
                    // if net rate is requested by user in chemDB, then add rate backward
                    if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                        
                        
                        topass = "prod";
                        double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                        
                        ri -= r0 * Qprod / chem->Keq[rxid] * exp((1.-ki)*(DUi)/KT ) ;
                        
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
                
                std::string msg;
                if (msk->wplog) {
                msg += "\n";
                output->toplog(msg);
                }
        
        DTtot*=(double)((nrL-1)*nrS);

        //=====================================================MICRO FAST RUN (END)==================================================================
        
        
        // Add to Dtot the contribution from the last layer
        //std::string msg;
        if (msk->wplog) {
            msg += "PARTICLE ";
            std::ostringstream ss;    ss << i;   msg = msg+ss.str()+"; Last layer; ";
        }
        
        //=============================================MICRO LAYER-BY-LAYER SLOW RUN (START)===========================================================
        /*for (int jj=0; jj<nrS; jj++){
            
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
                
                double Vt = (chem -> rx_dVt_fgd[rxid]);  // Tributary volume of fgd created by the reaction at the current step in the chain.
                
                // lattice cell volume: needed for rate
                double DV =  fix->fKMC_DV[pos];
                
                double DUi = 0.;
                if (strcmp(chem->mechinter[mid].c_str(),"int_no")!=0) {
                    DUi = DUpu;
                    if (chem->mechchain[mid]) DUi *= (chem->ch_rdV_fgd[chID][k]);
                    
                    // Add change in strain energy internal to the molecules
                    for (int j1=0; j1<chem->fgd_molID[rxid].size(); j1++){
                        int molid = chem->fgd_molID[rxid][j1];
                        DUi += chem->mol_Us[molid] * chem->mol_vm[molid] * chem->fgd_nmol[rxid][j1];
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
                
                r0 = r0*pow(Vt,(dim/3. - 1.))*DV;
                
                double ki = (chem -> ki[rxid]);
                
                // Forward rate equation: NB this includes contribution from largest interfacial energy change at contact with other phases (from function computing coverage areas)
                double uk_area;    // area of a unit kink (assuming cubic units)
                if (chem->mechchain[mid]) {  //if reaction is a chain
                    uk_area = pow(chem->ch_dVp_fgd[chID],2./3.) * 3.;
                }
                else{ //if instead it is a single reaction
                    uk_area = pow(chem->rx_dVp_fgd[rxid],2./3.) * 3.;
                }
                
                ri = r0 * Qreac * exp((ki)*( - DUi - tGM[i] * uk_area )/KT ) ;
                
                
                // if net rate is requested by user in chemDB, then add rate backward
                if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                    
                    topass = "prod";
                    double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    ri -= r0 * Qprod / chem->Keq[rxid] * exp((1.-ki) * ( DUi + tGM[i] * uk_area)/ KT );
                    
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
        }*/
        //=============================================MICRO LAYER-BY-LAYER SLOW RUN (END)===========================================================
        
        //=====================================================MICRO FAST RUN (START)==================================================================
        
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
                
                double Vt = (chem -> rx_dVt_fgd[rxid]);  // Tributary volume of fgd created by the reaction at the current step in the chain.
                
                // lattice cell volume: needed for rate
                double DV =  fix->fKMC_DV[pos];
                
                double DUi = 0.;
                if (strcmp(chem->mechinter[mid].c_str(),"int_no")!=0) {
                    DUi = DUpu;
                    if (chem->mechchain[mid]) DUi *= (chem->ch_rdV_fgd[chID][k]);
                    
                    // Add change in strain energy internal to the molecules
                    for (int j1=0; j1<chem->fgd_molID[rxid].size(); j1++){
                        int molid = chem->fgd_molID[rxid][j1];
                        DUi += chem->mol_Us[molid] * chem->mol_vm[molid] * chem->fgd_nmol[rxid][j1];
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
                
                r0 = r0*pow(Vt,(dim/3. - 1.))*DV;
                
                double ki = (chem -> ki[rxid]);
                
                // Forward rate equation: NB this includes contribution from largest interfacial energy change at contact with other phases (from function computing coverage areas)
                double uk_area;    // area of a unit kink (assuming cubic units)
                if (chem->mechchain[mid]) {  //if reaction is a chain
                    uk_area = pow(chem->ch_dVp_fgd[chID],2./3.) * 3.;
                }
                else{ //if instead it is a single reaction
                    uk_area = pow(chem->rx_dVp_fgd[rxid],2./3.) * 3.;
                }
                
                ri = r0 * Qreac * exp((ki)*( - DUi - tGM[i] * uk_area )/KT ) ;
                
                
                // if net rate is requested by user in chemDB, then add rate backward
                if (strcmp((chem->mechmode[mid]).c_str(),"net")==0) {
                    
                    topass = "prod";
                    double Qprod = solution->compQ(rxid,topass,fix->fKMCsinST[pos],fix->fKMCsinUL[pos]);
                    ri -= r0 * Qprod / chem->Keq[rxid] * exp((1.-ki) * ( DUi + tGM[i] * uk_area)/ KT );
                    
                    if (msk->wplog) {
                        std::ostringstream ss;
                        msg += "; ri ";    ss << ri;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Qprod ";    ss << Qprod;   msg += ss.str();    ss.str(""); ss.clear();
                        msg += "; Keq ";    ss << chem->Keq[rxid] ;   msg += ss.str();    ss.str(""); ss.clear();
                    }
                    
                }
                
                if (ri < 0.) ri = 0.;
                
                DTi = 1./ri;
                DTtot += DTi*(double)nrS;
                
                if (msk->wplog) {
                    msg += ", DT ";
                    std::ostringstream ss;    ss << DTi;   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", DTtot ";
                    ss << DTtot;   msg += ss.str();
                }
            }
            
            if (msk->wplog){
            msg+="\n";
            output->toplog(msg);
            }

        //=====================================================MICRO FAST RUN (END)==================================================================
        
        if (flag_bulk) rate_each[i] = 0.;
        else rate_each[i] = 1./DTtot;
        
        if (rate_each[i]<0.) rate_each[i] = 0.;    // if the backward ri's are > forward ri's the overall rate may end up < 0, which means that the current deletion event should not happen, hence its rate should be zero (not negative..)   CHECK THAT THIS DOES NOT GIVE PROBLEMS WHEN SELECTING THE EVENT TO CARRY OUT FROM CUMULATIVE RATE VECTORS, IN PARTICULAR WITH THE BINARY SEARCH ALGORITHM
        
        
        if (msk->wplog) {
            std::ostringstream ss;
            std::string msg = "\nPARTICLE ";
            ss << tID[i];   msg += ss.str();  ss.str("");   ss.clear();
            msg += " has nucleation rate  ";
            ss << rate_each[i];   msg += ss.str();
            output->toplog(msg);
        }
    }
    
}

/*
 
// ---------------------------------------------------------------
//  each processor computes its deletion rates
void Fix_delete::comp_rates_allser(int pos)
{
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
    
}



*/

// ---------------------------------------------------------------
//  the submaster receives rates from all slave processors and orders them in the rates vector corresponding to IDsrt (which is the only one that will be remembered)
void Fix_nucleate::rates_to_submaster(int pos){
    
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
double Fix_nucleate::CrateAll()
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
double Fix_nucleate::CratesTupdate(double Dt)
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
void Fix_nucleate::findEV(double CRC)
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
void Fix_nucleate::extractPROPS(int EVpID, int EVafixID)
{
    // this is called by all processors in subcomm..
   
    // shape and type parameters first, can be read directly from fix_nucleate command as they are user-given
    if (strcmp((fix->afKMCpgeom[EVafixID]).c_str(),"sphere")==0) {
        EVpDIAM =  fix->afKMCpdiam[EVafixID];
    }
    EVpTYPE =  fix->afKMCptype[EVafixID];

    
    // if ellipsoid instead... to be implemented
    
    // then position, to be read from lammps
    std::string tolmp;
    std::ostringstream typ;
    typ << EVpID;
    tolmp = "variable xptemp equal x["+typ.str()+"]";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable yptemp equal y["+typ.str()+"]";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable zptemp equal z["+typ.str()+"]";
    lammpsIO->lammpsdo(tolmp);
   
    tolmp = "run 0";
    lammpsIO->lammpsdo(tolmp);
    
    EVpXpos = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"xptemp",0));
    EVpYpos = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"yptemp",0));
    EVpZpos = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"zptemp",0));
   
    lammpsIO->lammpsdo("variable xptemp delete");
    lammpsIO->lammpsdo("variable yptemp delete");
    lammpsIO->lammpsdo("variable zptemp delete");
    
}





// ---------------------------------------------------------------
// Deleting selected particle with pID
void Fix_nucleate::execute(int pID, int EVafixID, int EVpTYPE,double EVpDIAM,double EVpXpos, double EVpYpos, double EVpZpos)
{
    if (msk->wplog) {
        std::string msg = "\n ****************** \n Nucleating particle ";
        std::ostringstream ss;
        ss << pID;   msg += ss.str()+"\n********************\n";
        output->toplog(msg);
    }
    
    
    // lammps deposits a particle at a random location within a tiny region around the particle target position (I am assuming that even a particle deposited out of the box will be ok because later moved inside. If not OK, in the future I might want to create region around centre of the simulation box, complex for triclinic ones though)
    std::string tolmp;
    std::ostringstream ss;
    
    double xl, xr, yl, yr, zl, zr;
    xl = EVpXpos - 1e-200;     xr = EVpXpos + 1e-200;
    yl = EVpYpos - 1e-200;     yr = EVpYpos + 1e-200;
    zl = EVpZpos - 1e-200;     zr = EVpZpos + 1e-200;
    
    /*
     tolmp = "region tempregg block ";
    ss << xl;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += " ";
    ss << xr;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += " ";
    ss << yl;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += " ";
    ss << yr;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += " ";
    ss << zl;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += " ";
    ss << zr;   tolmp += ss.str(); ss.str(""); ss.clear();
     */
    
    tolmp = "variable xltemp equal xlo+0.1*(xhi-xlo)";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable xrtemp equal xhi-0.1*(xhi-xlo)";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable yltemp equal ylo+0.1*(yhi-ylo)";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable yrtemp equal yhi-0.1*(yhi-ylo)";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable zltemp equal zlo+0.1*(zhi-zlo)";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable zrtemp equal zhi-0.1*(zhi-zlo)";
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "region tempregg block $(v_xltemp) $(v_xrtemp) $(v_yltemp) $(v_yrtemp) $(v_zltemp) $(v_zrtemp) units box";
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "group gtempin empty";
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "fix fix_dep_temp gtempin deposit 1 ";
    ss << EVpTYPE;   tolmp += ss.str(); ss.str(""); ss.clear();
    tolmp += " 1 12345 region tempregg";    //by default, max id is assigned, and particle is accepted unless exactly on top of another existing one
    lammpsIO->lammpsdo(tolmp);
    
    lammpsIO->lammpsdo("run 1");
    
    tolmp = "variable xltemp delete";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable xrtemp delete";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable yltemp delete";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable yrtemp delete";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable zltemp delete";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "variable zrtemp delete";
    lammpsIO->lammpsdo(tolmp);
    
    // read id of new particle, and assign it the exact desired position
    tolmp = "compute cidtemp gtempin property/atom id";
    lammpsIO->lammpsdo(tolmp);
    //lammpsIO->lammpsdo("run 1");
    
    //temp dump to update variables and computes
    tolmp = "dump tdID all custom 1 dump.tempIN_"+universe->SCnames[universe->color]+" id type x y z radius c_cidtemp";
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "compute cidtt gtempin reduce sum c_cidtemp";
    lammpsIO->lammpsdo(tolmp);
    
    lammpsIO->lammpsdo("run 1");
    lammpsIO->lammpsdo("unfix fix_dep_temp");
    
    //natoms = static_cast<int> (lammpsIO->lmp->atom->natoms); // total number of atoms in lammps, all groups all fixes)
    //tolmp = "variable natt equal ;
    //lammpsIO->lammpsdo(tolmp);
    
    
    // alters lammps thermo to compute the required quantities
    tolmp = lammpsIO->lmpThSt + " c_cidtt";
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");
    
    
    //tolmp = "variable vidtt equal $(c_cidtemp[";
    //ss << natoms;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += "])";
    tolmp = "variable vidtt equal $(c_cidtt)";
    lammpsIO->lammpsdo(tolmp);

   
    // alters lammps thermo to compute the required quantities
    tolmp = lammpsIO->lmpThSt + " c_cidtt v_vidtt";
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");
    
    
    //tolmp = "variable vidtt equal $(c_cidtemp[1])";
    //tolmp = "variable vidtt equal $(c_cidtt)";
    //lammpsIO->lammpsdo(tolmp);

    tolmp = "set atom $(v_vidtt) diameter ";
    ss << EVpDIAM;   tolmp += ss.str(); ss.str(""); ss.clear();
    //fprintf(screen,"\n PROC %d : \"%s\"\n",me,tolmp.c_str());
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 0");

    // fprintf(screen,"DEBUG 1: PROC %d : \"%s\"\n",me,tolmp.c_str());
    // sleep (1);
    
    tolmp = "set atom $(v_vidtt) x ";
    ss << EVpXpos;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += " y ";
    ss << EVpYpos;   tolmp += ss.str(); ss.str(""); ss.clear(); tolmp += " z ";
    ss << EVpZpos;   tolmp += ss.str(); ss.str(""); ss.clear();
    //fprintf(screen,"\n PROC %d : \"%s\"\n",me,tolmp.c_str());
    lammpsIO->lammpsdo(tolmp);


    // std::ostringstream typ;
    // typ << fix->fKMCptypeTRY[pos];
    // tolmp = "group gtemp type "+typ.str();
    // //tolmp = "group gtemp type 4";
    // lammpsIO->lammpsdo(tolmp);
    // tolmp = "displace_atoms gtemp random 1E-8 1E-8 1E-8 12345 units box";
    // lammpsIO->lammpsdo(tolmp);

    tolmp = "group gtemp type "; 
    ss << EVpTYPE;   tolmp += ss.str(); ss.str(""); ss.clear();
    lammpsIO->lammpsdo(tolmp);
    tolmp = "displace_atoms gtemp random 1E-20 1E-20 1E-20 12345 units box";
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("group gtemp delete");

    // fprintf(screen,"DEBUG 2: PROC %d : \"%s\"\n",me,tolmp.c_str());
    // sleep (1);    
    lammpsIO->lammpsdo("run 0");

    // fprintf(screen,"DEBUG 3: PROC %d : \"%s\"\n",me,tolmp.c_str());
    // sleep (1);
    
    // tolmp = "set atom $(v_vidtt) diameter ";
    // ss << EVpDIAM;   tolmp += ss.str(); ss.str(""); ss.clear();
    // //fprintf(screen,"\n PROC %d : \"%s\"\n",me,tolmp.c_str());
    // lammpsIO->lammpsdo(tolmp);
    // lammpsIO->lammpsdo("run 0");

    

    tolmp = "undump tdID";     // removing temporary dump
    lammpsIO->lammpsdo(tolmp);
    
    
    //lammpsIO->lammpsdo("unfix fix_dep_temp");
    lammpsIO->lammpsdo("variable vidtt delete");
    
    lammpsIO->lammpsdo(lammpsIO->lmpThSt);    //restore original lammps thermo
    
    lammpsIO->lammpsdo("uncompute cidtt");
    lammpsIO->lammpsdo("uncompute cidtemp");
    lammpsIO->lammpsdo("region tempregg delete");
    lammpsIO->lammpsdo("group gtempin clear");
    lammpsIO->lammpsdo("group gtempin delete");

   
    //fprintf(screen,"\n PROC %d restoring lammps thermo command \"%s\"\n",me,lammpsIO->lmpThSt.c_str());
    //   sleep(2);
    
    //lammpsIO->lammpsdo("variable natt delete");
    
    
    if (msk->wplog) {
        std::string msg = "\n Succesfully created particle in lammps \n";
        output->toplog(msg);
    }
            
    
    //tolmp = "compute tempRAD sub property/atom radius";
    //lammpsIO->lammpsdo(tolmp);
    //double *tempRad;
    //tempRad = new double[1];
    //tempRad = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempRAD",1,1));
    
    //tempRad = ((double *) lammps_extract_atom(lammpsIO->lmp,(char *)"radius"));
    
    
    partVol =  M_PI / 6. * EVpDIAM*EVpDIAM*EVpDIAM;  // particle volume  FOR SPHERICAL
    
    //tolmp = "uncompute tempRAD";
    //lammpsIO->lammpsdo(tolmp);
    
    /*
     tolmp = "delete_atoms group sub";
    lammpsIO->lammpsdo(tolmp);
    tolmp = "group sub delete";
    lammpsIO->lammpsdo(tolmp);
     */
    //tolmp = "run 1";   // just to see something in dump... later to be replaced with invoking a fix_relax
    //lammpsIO->lammpsdo(tolmp);
    
}






// ---------------------------------------------------------------
//  All procs (although only submaster is important...) adds elements to cumulative vector of number of particles in fix, and associated absolute ID of fix
void Fix_nucleate::add_Ng_apos(int Np, int apos)
{
    if (cumNPfix.size()==0)cumNPfix.push_back(Np);
    else cumNPfix.push_back(cumNPfix.back()+Np);
    
    fixaID.push_back(apos);
}






// ---------------------------------------------------------------
// Each submaster builds vectors with particle IDs, energies, rates etc, from all procs in subcomm. This function clears them
void Fix_nucleate::clearvecs()
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
void Fix_nucleate::printall()
{
	fprintf(screen,"\n---------ALL ABOUT NUCLEATE----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
