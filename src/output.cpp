#include "output.h"
#include "chemistry.h"
#include "universe.h"
#include "lammpsIO.h"
#include "error.h"

#include <string.h>

using namespace MASKE_NS;

// ---------------------------------------------------------------
// Initialize class
Output::Output(MASKE *maske) : Pointers(maske)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating output class\n");
    
    th_every = 0;
    
}




// ---------------------------------------------------------------
// Class destructor
Output::~Output()
{
    
}



// ---------------------------------------------------------------
// creates processor specific log file
void Output::createplog(std::string fname)
{
    plog = fopen (fname.c_str(),"w");
    fclose(plog);
    toplog("BEGINNING OF PROCESSOR SPECIFIC LOG");
    
}



// ---------------------------------------------------------------
// Write to processor-specific log file
void Output::toplog(std::string msg)
{
    plog = fopen (msk->plogfname.c_str(),"a");
    fprintf(plog,"%s\n",msg.c_str());
    fclose(plog);
}



// ---------------------------------------------------------------
// create thermo file (master only)
void Output::createthermo(std::string fname)
{
    thermo = fopen (fname.c_str(),"w");
    fprintf(thermo,"Step\tTime");
    for (int i=0; i<th_qtts.size(); i++) fprintf(thermo,"\t%s",th_qtts[i].c_str());
    fprintf(thermo,"\n");
    fclose(thermo);
    
}


// ---------------------------------------------------------------
// add quantity to output in thermo file
void Output::add_thqtt(std::string name)
{
    th_qtts.push_back(name);
}


// ---------------------------------------------------------------
// write new entry in thermo file (master only)
void Output::writethermo(void)
{
    if (me == MASTER){
        thermo = fopen (msk->th_fname.c_str(),"a");
        fprintf(thermo,"%d",msk->step);
        fprintf(thermo,"\t%e",msk->tempo);
    }
    
    for (int i=0; i<th_qtts.size(); i++) {
        std::istringstream iss(th_qtts[i]);
        std::string token;
        std::getline(iss, token, '_');   // reads whatever before the first underscore
        if (strcmp(token.c_str(),"conc")==0) {
            iss>>token; // now this contains molecule name
            // locate molecule position in chem
            for (int j=0; j<chem->Nmol; j++) {
                if (strcmp((chem->molnames[j]).c_str(),token.c_str())==0) {
                    if (me == MASTER) fprintf(thermo,"\t%e",chem->mol_cins[j]);
                }
            }
        } else if (strcmp(token.c_str(),"lmp")==0) {
            std::string vc; //either v or c, to identify lammps variable or compute
            std::getline(iss, vc, '_');
            iss>>token; // now this contains a variable or compute name from lammps
            double varc = 0.;
            
            // the first subcomm with lammps active records the lmp variable or compute value
            if (universe->flampSC == universe->color){
                std::string tolmp = "run 0";
                lammpsIO->lammpsdo(tolmp);
                 if (strcmp(vc.c_str(),"v")==0) {
                     fprintf(screen,"PROC %d varname %s",me,token.c_str());
                     varc = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)token.c_str(),0));
                 }
                 else if (strcmp(vc.c_str(),"c")==0) {
                     varc = *((double *) lammps_extract_compute(lammpsIO->lmp,(char *)token.c_str(),0,0));
                 }
            }
        
            // If first subcomm with lammps active is not the MASTER's one, communicate value to master
            if (universe->flampSC != 0){
            }
            //      if universe->flampSC == universe->color   && key==0
            //          send variable value out
            //      else
            //          recv from universe->flampID
            
            if (me == MASTER) fprintf(thermo,"\t%e",varc);
        }
    }
    if (me == MASTER) fprintf(thermo,"\n");
    if (me == MASTER) fclose(thermo);
    
}



// ---------------------------------------------------------------
// record new dump
void Output::add_dump(std::string dID,int devery,std::string dstring)
{
    dumpID.push_back(dID);
    dump_every.push_back(devery);
    dstring = "dump "+dID+" "+dstring;
    dump_string.push_back(dstring);
    dump_first.push_back(true);
   
    
    
    // find numeric ID of tempdump associated to this dump
    /*int tdIDnum = -2;
    for (int i =0; i<tempdumpID.size(); i++) {
        if (strcmp(tdID.c_str(),tempdumpID[i].c_str())==0) {
            tdIDnum = i;
        }
    }
    if (strcasecmp(tdID.c_str(),"null")==0) {
        tdIDnum = -1;
    }
    if (tdIDnum<-1) {
        std::string msg = "ERROR: tempdump "+tdID+" invoked by dump \""+dID+"\" not found.";
        error->errsimple(msg);
    }
    else {
        TdumpID.push_back(tdIDnum);
    }
     */
}


// ---------------------------------------------------------------
// write new entry in dump file
void Output::writedump(int i)
{
    /*if (TdumpID[i]!=-1) {
        openTdump(TdumpID[i]);
        lammpsIO->lammpsdo("run 1");
    }*/
    
    
    // create dump (for just one step)
    lammpsIO->lammpsdo(dump_string[i]);
    
    // make it run immediately and append if not first
    std::string tolmp;
    if (dump_first[i]) {
        tolmp = "dump_modify "+dumpID[i]+" every 1 first yes";
        lammpsIO->lammpsdo(tolmp);
    }
    else {
        tolmp = "dump_modify "+dumpID[i]+" every 1 first yes append yes";
        lammpsIO->lammpsdo(tolmp);
    }
    
    // write enrty
    lammpsIO->lammpsdo("run 0");
    
    // close dump
    tolmp = "undump "+dumpID[i];
    lammpsIO->lammpsdo(tolmp);
    
    //if (TdumpID[i]!=-1) closeTdump(TdumpID[i]);
}





// ---------------------------------------------------------------
// record new temporary dump
void Output::add_tempdump(std::string dID,std::string dstring)
{
    tempdumpID.push_back(dID);
    tempdump_string.push_back(dstring);
    std::string modif;
    modif = "dump_modify "+dID+" first yes every 1";
    tempdump_modify.push_back(modif);
}


// ---------------------------------------------------------------
// sets up a tempdump
void Output::openTdump(int i)
{
    lammpsIO->lammpsdo(tempdump_string[i]);
    lammpsIO->lammpsdo(tempdump_modify[i]);
}


// ---------------------------------------------------------------
// closes a tempdump
void Output::closeTdump(int i)
{
    std::string tolmp;
    tolmp = "undump "+tempdumpID[i];
    lammpsIO->lammpsdo(tolmp);
}

