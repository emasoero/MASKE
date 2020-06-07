#include "fix.h"
#include <sstream>
#include "universe.h"
#include "error.h"
#include "chemistry.h"
#include "block.h"
#include "store.h"
//#include "lammpsIO.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#include <string.h>

using namespace MASKE_NS;

Fix::Fix(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    //Ntrans = 0;
    //Ntsc=0;

}



// ---------------------------------------------------------------
// Class destructor
Fix::~Fix()
{
    
}



// ---------------------------------------------------------------
// Add fix to list of active fixes
void Fix::add(std::string instr)
{
    std::string PROCtype;  //process type can be fKMC, Krelax, or Cont
    std::istringstream ss(instr);
    
    std::string temp_type, temp_name, temp_scom,keyword, temp_mech,temp_wtype, temp_sinST , temp_sinUL, temp_sinbox , temp_soutUL ,temp_soutbox, temp_reg, temp_latt, temp_pgeom, temp_min, temp_store, temp_group;
    int temp_ptype, temp_ptypeTRY;
    double temp_pdiam;

    //some default values follow
    temp_sinST = "fixed";
    temp_sinUL = "uniform";
    temp_sinbox = "box";
    temp_soutUL = "uniform";
    temp_soutbox = "box";
    double temp_evt=0., temp_leval=0., temp_dt=0.;
    int temp_mid=-1, temp_rid=-1, temp_lid=-1, temp_minid=-1, temp_sid=-1, temp_steps=-1;

    
    std::vector<double> wargs;
    
    
    
    // all processors read the fix string but later only the MASTER and the invoked subcomm keep track of them
    ss >> PROCtype;
    if (strcmp(PROCtype.c_str(),"KMC-free")==0) {
        ss >> temp_type;
        ss >> temp_name;
        check_name_unique(temp_name);
        allFixNames.push_back(temp_name);
        
        ss >> temp_scom;
        if (strcmp(temp_type.c_str(),"delete")==0) {
            ss >> temp_ptype;  // the only useful parameter for delete: all that follows are defaults
            temp_reg = "none";
            temp_rid = -1;
            temp_latt = "none";
            temp_lid = -1;
            temp_min = "none";
            temp_minid = -1;
            temp_ptypeTRY = -1;
            temp_pgeom = "none";
            temp_pdiam = -1;
        }
        else if (strcmp(temp_type.c_str(),"nucleate")==0) {
            ss >> temp_reg >> temp_latt >> temp_min >> temp_ptypeTRY >> temp_ptype >> temp_pgeom;
            
            for (int i=0; i<store->Nreg; i++) {
                if (strcmp(temp_reg.c_str(),(store->RegNames[i]).c_str())==0) temp_rid=i;
            }
            if (temp_rid==-1) {
                std::string msg = "ERROR: region name \""+temp_reg+"\" specified for fix \""+temp_name+"\" does not match any stored region. \n";
                error->errsimple(msg);
            }
            
            for (int i=0; i<store->Nlat; i++) {
                if (strcmp(temp_latt.c_str(),(store->LatNames[i]).c_str())==0) temp_lid=i;
            }
            if (temp_lid==-1) {
                std::string msg = "ERROR: lattice name \""+temp_latt+"\" specified for fix \""+temp_name+"\" does not match any stored lattice. \n";
                error->errsimple(msg);
            }
            
            for (int i=0; i<store->Nmin; i++) {
                if (strcmp(temp_min.c_str(),(store->MinNames[i]).c_str())==0) temp_minid=i;
            }
            if (temp_minid==-1) {
                std::string msg = "ERROR: minimiser name \""+temp_min+"\" specified for fix \""+temp_name+"\" does not match any stored minimiser. \n";
                error->errsimple(msg);
            }
            
            if (strcmp(temp_pgeom.c_str(),"sphere")==0) ss >> temp_pdiam;
            else {
                std::string msg = "ERROR: fix \""+temp_name+"\" is trying to nucleate a particle with unknown geometry \""+temp_pgeom+"\". Only \"sphere\" is supported as a geometry to date.\n";
                error->errsimple(msg);
            }
        }
        else {
            std::string msg = "ERROR: fix \""+temp_name+"\" is of unknown type \""+temp_type+"\". Only \"delete\" and \"nucleate\" are supported types to date.\n";
            error->errsimple(msg);
        }
        
        bool commt_found = false;
        while (!ss.eof() && !commt_found) {
            ss >> keyword;
            if (strcmp(keyword.c_str(),"mech")==0) {
                ss >> temp_mech;
                // search for mech in chemistry database
                for (int i=0; i<chem->Nmech; i++) {
                    if (strcmp(temp_mech.c_str(),(chem->mechnames[i]).c_str())==0) temp_mid=i;
                }
                if (temp_mid==-1) {
                    std::string msg = "ERROR: mechanism name \""+temp_mech+"\" specified for fix \""+temp_name+"\" does not match with any of the mechanism names specificed in the chemistry data base. \n";
                    error->errsimple(msg);
                }
            }
            else if (strcmp(keyword.c_str(),"everyt")==0) {
                ss >> temp_evt;
            }
            else if (strcmp(keyword.c_str(),"leval")==0) {
                ss >> temp_leval;
            }
            else if (strcmp(keyword.c_str(),"wei")==0) {
                ss >> temp_wtype;
                int nargs = 0;
                if (strcmp(temp_wtype.c_str(),"simple")==0) nargs = 1;
                else{
                    std::string msg = "ERROR: invalid weight calculatore type in fix "+instr+" \n \""+temp_wtype+"\" not recognized \n";
                    error->errsimple(msg);
                }
                for(int i=0; i<nargs; i++) {
                    double targ;
                    ss >> targ;
                    wargs.push_back(targ);
                }
            }
            else if (strcmp(keyword.c_str(),"sol_in")==0) {
                ss >> temp_sinST >> temp_sinUL >> temp_sinbox ;
            }
            else if (strcmp(keyword.c_str(),"sol_out")==0) {
                ss >> temp_soutUL >> temp_soutbox;
            }
            else if (keyword[0]=='#'){
                commt_found = true;
            }
            else {
                std::string msg = "ERROR: invalid keyword in fix "+temp_name+" \n "+keyword+" not recognized \n";
                error->errsimple(msg);
            }
        }
        
        if (wargs.size()==0) {
            temp_wtype = "simple";
            wargs.push_back(1.);    //default weight calculator
        }
        
        
        // Adding fix to all-fixes list: only the MASTER will use it in krun, but I let all processors know it for when they will have to locate each local event in the global list: this will reduce communication
        afKMCtype.push_back(temp_type);
        afKMCname.push_back(temp_name);
        afKMCscom.push_back(temp_scom);
        
        afKMCreg.push_back(temp_reg);
        afKMCrid.push_back(temp_rid);
        afKMClatt.push_back(temp_latt);
        afKMClid.push_back(temp_lid);
        afKMCmin.push_back(temp_min);
        afKMCminid.push_back(temp_minid);
        afKMCptypeTRY.push_back(temp_ptypeTRY);
        afKMCptype.push_back(temp_ptype);
        afKMCpgeom.push_back(temp_pgeom);
        afKMCpdiam.push_back(temp_pdiam);
        
        afKMC_DV.push_back(0.);
        
        afKMCmech.push_back(temp_mech);
        afKMCmid.push_back(temp_mid);
        afKMCeveryt.push_back(temp_evt);
        afKMCleval.push_back(temp_leval);
        afKMCnevents.push_back(0);  //to be computed dynamically during krun
        afKMCcumRate.push_back(0.);  // to be computed by each fix_ invoked by krun
        afKMCwei.push_back(temp_wtype);
        afKMCwarg.push_back(wargs);
        afKMCsinST.push_back(temp_sinST);
        afKMCsinUL.push_back(temp_sinUL);
        afKMCsinbox.push_back(temp_sinbox);
        afKMCsoutUL.push_back(temp_soutUL);
        afKMCsoutbox.push_back(temp_soutbox);
        
        
    
        bool added_local_fix = false;
        //inndividual processors in fix-invoked subcomm adding fix to their local list
        if (strcmp(temp_scom.c_str(),(universe->SCnames[universe->color]).c_str())==0) {
            fKMCaID.push_back(afKMCname.size()-1);
            fKMCtype.push_back(temp_type);
            fKMCname.push_back(temp_name);
            fKMCscom.push_back(temp_scom);
            
            fKMCreg.push_back(temp_reg);
            fKMClatt.push_back(temp_latt);
            fKMCptypeTRY.push_back(temp_ptypeTRY);
            fKMCptype.push_back(temp_ptype);
            fKMCpgeom.push_back(temp_pgeom);
            fKMCpdiam.push_back(temp_pdiam);
            fKMCmin.push_back(temp_min);
            fKMCminid.push_back(temp_minid);
            
            fKMC_DV.push_back(0.);
            
            fKMCmech.push_back(temp_mech);
            fKMCmid.push_back(temp_mid);
            fKMCrid.push_back(temp_rid);
            fKMClid.push_back(temp_lid);
            fKMCeveryt.push_back(temp_evt);
            fKMCleval.push_back(temp_leval);
            fKMCnevents.push_back(0);
            fKMCwei.push_back(temp_wtype);
            fKMCwarg.push_back(wargs);
            fKMCsinST.push_back(temp_sinST);
            fKMCsinUL.push_back(temp_sinUL);
            fKMCsinbox.push_back(temp_sinbox);
            fKMCsoutUL.push_back(temp_soutUL);
            fKMCsoutbox.push_back(temp_soutbox);
            fKMCglobID.push_back(afKMCtype.size()-1);
            fKMCrank.push_back(-1);  // to be computed in krun
            fKMCfirst.push_back(-1);  // to be computed in krun
            added_local_fix = true;
        }
        
        // Each submaster prints their recorded fix, just to check that all is fine
        if (universe->key==0 && added_local_fix) {
            int siz;
            siz = (int)fKMCtype.size()-1;
            
            //fprintf(screen,"\n mech id is %d ",fKMCmid[siz]);
            
            if (strcmp(fKMCtype[siz].c_str(),"delete")==0)
            {
                fprintf(screen,"\n Added fix: KMC-free %s %s %s %d mech %s everyt %f sol_in_style %s sol_in_unif %s sol_in_box %s sol_out_unif %s sol_out_box %s wei %s \n", fKMCtype[siz].c_str(),fKMCname[siz].c_str(),fKMCscom[siz].c_str(),fKMCptype[siz],(chem->mechnames[fKMCmid[siz]]).c_str(),fKMCeveryt[siz],fKMCsinST[siz].c_str(),fKMCsinUL[siz].c_str(),fKMCsinbox[siz].c_str(),fKMCsoutUL[siz].c_str(),fKMCsoutbox[siz].c_str(),fKMCwei[siz].c_str());
            }
            else if (strcmp(fKMCtype[siz].c_str(),"nucleate")==0){
                    fprintf(screen,"\n Added fix: KMC-free %s %s %s %s %s %s %d %d %s %f mech %s everyt %f sol_in_style %s sol_in_unif %s sol_in_box %s sol_out_unif %s sol_out_box %s wei %s \n", fKMCtype[siz].c_str(),fKMCname[siz].c_str(),fKMCscom[siz].c_str(),(store->RegNames[fKMCrid[siz]]).c_str(),(store->LatNames[fKMClid[siz]]).c_str(),(store->MinNames[fKMClid[siz]]).c_str(),fKMCptypeTRY[siz],fKMCptype[siz],fKMCpgeom[siz].c_str(),fKMCpdiam[siz],(chem->mechnames[fKMCmid[siz]]).c_str(),fKMCeveryt[siz],fKMCsinST[siz].c_str(),fKMCsinUL[siz].c_str(),fKMCsinbox[siz].c_str(),fKMCsoutUL[siz].c_str(),fKMCsoutbox[siz].c_str(),fKMCwei[siz].c_str());
            }
            int nargs = 0;
            if (strcmp(temp_wtype.c_str(),"simple")==0) nargs = 1;
            for(int i=0; i<nargs; i++) fprintf(screen,"%f ", fKMCwarg[siz][i]);
            fprintf(screen,"\n Local fix position: %d    Global fix position: %d \n",siz,fKMCglobID[siz]);
        }
    }
    
    
    else if (strcmp(PROCtype.c_str(),"KMC-rej")==0) {
        ss >> temp_type;
        ss >> temp_name;
        check_name_unique(temp_name);
        allFixNames.push_back(temp_name);

    }
    else if (strcmp(PROCtype.c_str(),"Cont")==0) {
        ss >> temp_type;
        ss >> temp_name;
        check_name_unique(temp_name);
        allFixNames.push_back(temp_name);
        ss >> temp_scom;
	if (temp_type == "nufeb") {
	  ss >> temp_store;
	  temp_sid = -1;
	  for (int i=0; i<store->MulNames.size(); i++) {
	    if (store->MulNames[i] == temp_store)
	      temp_sid = i;
	  }
	  if (temp_sid < 0) {
	    error->errsimple("ERROR: couldn't find multiple-line store for nufeb fix");
	  }
	} else {
	  temp_store = "none";
	}
        bool commt_found = false;
	temp_steps = -1;
        while (!ss.eof() && !commt_found) {
            ss >> keyword;
            if (strcmp(keyword.c_str(),"leval")==0) {
                ss >> temp_leval;
            }
            else if (strcmp(keyword.c_str(),"dt")==0) {
                ss >> temp_dt;
            }
            else if (strcmp(keyword.c_str(),"steps")==0) {
                ss >> temp_steps;
            }
            else if (strcmp(keyword.c_str(),"group")==0) {
                ss >> temp_group;
            }
            else if (strcmp(keyword.c_str(),"sol_out")==0) {
                ss >> temp_soutbox;
            }
            else if (keyword[0]=='#'){
                commt_found = true;
            }
            else {
                std::string msg = "ERROR: invalid keyword in fix "+temp_name+" \n "+temp_wtype+" not recognized \n";
                error->errsimple(msg);
            }
        }
        
        // Adding fix to all-Cont-fixes list
        aCtype.push_back(temp_type);
        aCname.push_back(temp_name);
        aCscom.push_back(temp_scom);
        aCdt.push_back(temp_dt);
	aCleval.push_back(temp_leval);
	aCstore.push_back(temp_store);
	aCsid.push_back(temp_sid);
	aCsteps.push_back(temp_steps);
	aCgroups.push_back(temp_group);
	aCsoutbox.push_back(temp_soutbox);

        bool added_local_fix = false;
        //inndividual processors in fix-invoked subcomm adding fix to their local list
        if (strcmp(temp_scom.c_str(),(universe->SCnames[universe->color]).c_str())==0) {
            Ctype.push_back(temp_type);
            Cname.push_back(temp_name);
            Cscom.push_back(temp_scom);
            Cdt.push_back(temp_dt);
            Cleval.push_back(temp_leval);
	    Cstore.push_back(temp_store);
	    Csid.push_back(temp_sid);
	    Csteps.push_back(temp_steps);
	    Cgroups.push_back(temp_group);
	    Csoutbox.push_back(temp_soutbox);
            CglobID.push_back(aCtype.size()-1);
            added_local_fix = true;
        }
        
        // Each submaster prints their recorded fix, just to check that all is fine
        if (universe->key==0 && added_local_fix) {
            int siz;
            siz = (int)Ctype.size()-1;
            
            fprintf(screen,"\n Added fix: Cont %s %s %s %f", Ctype[siz].c_str(),Cname[siz].c_str(),Cscom[siz].c_str(),Cdt[siz]);
            int nargs = 0;
            fprintf(screen,"\n Local position: %d    Global position: %d",siz,CglobID[siz]);
        }
      

    }
    else if (strcmp(PROCtype.c_str(),"Krelax")==0) {
        ss >> temp_type;
        ss >> temp_name;
        check_name_unique(temp_name);
        allFixNames.push_back(temp_name);
    }
    else {
        std::string msg = "ERROR: Unknown process type. Valid process types are KMC-free, KMC-rej, Krelax, and Cont. Instead found "+ PROCtype+" \n";
        error->errsimple(msg);
    }
        
    return;
    //fprintf(screen,"\n String passed to fix_add is: \n %s \n",instr.c_str());
    //exit(0);
}




// ---------------------------------------------------------------
// Check that the fix name is unique. If not, error and exit
void Fix::check_name_unique(std::string temp_name)
{
    if (allFixNames.size()>0){
        for (int i=0; i<allFixNames.size(); i++) {
            if (strcmp(allFixNames[i].c_str(),temp_name.c_str())==0) {
                std::ostringstream id1, id2;
                id1 << i;
                id2 << allFixNames.size();
                std::string msg = "ERROR: global fixes number "+id1.str()+" and "+id2.str()+" have same name: "+temp_name+" \n";
                error->errsimple(msg);
            }
        }
    }
    
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Fix::printall()
{
	fprintf(screen,"\n---------ALL ABOUT FIX----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
