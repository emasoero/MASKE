#include "chemistry.h"
#include "memory.h"
#include "error.h"


//#include "universe.h"
//#include "DTnucleate.h"
//
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "clinker.h"
#include "pyramids.h"
#include "hydration.h"
#include "csh.h"
#include "c2s.h"
#include "c3s.h"
#include "water.h"
*/
#include <string.h>

using namespace MASKE_NS;


Chemistry::Chemistry(MASKE *maske) : Pointers(maske)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating molecules and reactions class\n");
    
    Nmol = 0;
    Ngx = 0;
    Nreax = 0;
    NDG = 0;
    Nreax = 0;
    Nchain = 0;
    Nsen = 0;
    Nmech=0;
}



// ---------------------------------------------------------------
// Class destructor
Chemistry::~Chemistry()
{
    //delete lmp;
}





// ---------------------------------------------------------------
// Read chemistry database file
void Chemistry::readDB(std::string fname)
{
    std::ifstream inCDB(fname.c_str());
    if (!inCDB.is_open()) {
        std::string msg = "ERROR: Chemistry database file  \""+fname+"\" not found.";
        error->errsimple(msg);
    }
    else {
        //int count_debug = 0;
        while (!inCDB.eof()) {
            std::string line, word;
            std::getline (inCDB, line);
            ss.clear();
            ss.str(line);
            ss >> word;
            if (strncmp(word.c_str(), "#", 1) == 0) int foo = 1; // Do nothing. Lines starting with # are comments
            else if (strcmp(word.c_str(), "molecule") == 0)   addmolecule();
            else if (strcmp(word.c_str(), "gammax") == 0)   addgammax();
            else if (strcmp(word.c_str(), "DGx") == 0)   addDG();
            else if (strcmp(word.c_str(), "reax") == 0)   addreax();
            else if (strcmp(word.c_str(), "sen") == 0)   addsen();
            else if (strcmp(word.c_str(), "mech") == 0)   addmech();
            else if (!word.empty()){
                std::string msg = "\nERROR: Unknown command in chemistry database file: "+word+"\n\n";
                error->errsimple(msg);
            }
            
            //count_debug  ++;
            //fprintf(screen," PROC %d Read Chem DB: %d \n",me, count_debug);
            
        }
    }
    
    
    
}




// ---------------------------------------------------------------
// Add molecule type as defined from input
void Chemistry::addmolecule()
{
    Nmol++;
    
    std::string newname;
    ss >> newname;
    molnames.push_back(newname);
    
    double size;
    ss >> size;
    mol_arad.push_back(size);
    ss >> size;
    mol_acir.push_back(size);
    ss >> size;
    mol_vapp.push_back(size);
    ss >> size;
    mol_ahyd.push_back(size);
    ss >> size;
    mol_bDH.push_back(size);
    ss >> size;
    mol_z.push_back(size);
    // We purposely include this even when MASKE_WITH_NUFEB is not defined because
    // we didn't want to change the input files if compiled with or without NUFEB
    int n;
    ss >> n;
    mol_nufeb.push_back(n);
    ss >> n;
    mol_nufeb_form.push_back(n);
    // Same as above for speciation
    std::string elem; 
    ss >> elem;
    mol_spec.push_back(elem);
    
    mol_cins.push_back(0.);  //initial concentration set to zero and later to be updated by solution.cpp
    mol_nins.push_back(0.);  //initial number of mols set to zero and later to be updated by solution.cpp
    mol_cindV.push_back(0.);  //initial concentration set to zero and later to be updated by solution.cpp
    mol_nindV.push_back(0.);  //initial number of mols set to zero and later to be updated by solution.cpp
    
    
      //just a check
    if (me==MASTER) fprintf(screen," New molecule added: %s %f %f %f %f %f\n",molnames[Nmol-1].c_str(),mol_arad[Nmol-1],mol_acir[Nmol-1],mol_vapp[Nmol-1],mol_ahyd[Nmol-1],mol_z[Nmol-1]);
}




// ---------------------------------------------------------------
// Add rules to calculate activity coefficient of activated complex
void Chemistry::addgammax()
{
    
     Ngx++;
    
    std::string newname;
    ss >> newname;
    gxnames.push_back(newname);
    
    ss >> newname;
    gxstyle.push_back(newname);
    
    int ncoef = 0;
    if (strcmp(newname.c_str(), "const") == 0) ncoef = 1;
    else if(strcmp(newname.c_str(), "test3coef") == 0) ncoef = 3;
    else {
        std::string msg = "ERROR: Unknown style of gammax calculator: \""+newname+"\" \n";
        error->errsimple(msg);
    }
    
    std::vector<double> cvec;
    for (int i=0; i<ncoef; i++) {
        double entry;
        ss >> entry;
        cvec.push_back(entry);
    }
    
    gxcoef.push_back(cvec);
    
    
    //just a check
    if (me==MASTER) fprintf(screen," New gammax calculator added: %s Style = %s\n",gxnames[Ngx-1].c_str(),gxstyle[Ngx-1].c_str());
    if (me==MASTER) fprintf(screen," Current gammax coeff matrix is: \n");
    for (int i=0; i<Ngx; i++) {
        for (int j=0; j<gxcoef[i].size(); j++) {
             if (me==MASTER) fprintf(screen," %f ",gxcoef[i][j]);
        }
         if (me==MASTER) fprintf(screen,"\n");
    }
}






// ---------------------------------------------------------------
// Add rules to calculate activity coefficient of activated complex
void Chemistry::addDG()
{
    
    
    NDG++;
    
    std::string newname;
    ss >> newname;
    DGnames.push_back(newname);
    
    ss >> newname;
    DGstyle.push_back(newname);
    
    int ncoef = 0;
    if (strcmp(newname.c_str(), "const") == 0) ncoef = 1;
    else if(strcmp(newname.c_str(), "test3coef") == 0) ncoef = 3;
    else {
        std::string msg = "ERROR: Unknown style of DGx calculator: \""+newname+"\" \n";
        error->errsimple(msg);
    }
    
    std::vector<double> cvec;
    for (int i=0; i<ncoef; i++) {
        double entry;
        ss >> entry;
        cvec.push_back(entry);
    }
    
    DGcoef.push_back(cvec);
    
    
    //just a check
    if (me==MASTER) fprintf(screen," New DGx calculator added: %s Style = %s\n",DGnames[NDG-1].c_str(),DGstyle[NDG-1].c_str());
    if (me==MASTER) fprintf(screen," Current DGx coeff matrix is: \n");
    for (int i=0; i<NDG; i++) {
        for (int j=0; j<DGcoef[i].size(); j++) {
            if (me==MASTER) fprintf(screen," %f ",DGcoef[i][j]);
        }
        if (me==MASTER) fprintf(screen,"\n");
    }
     
    
}





// ---------------------------------------------------------------
// Add reaction as defined from input
void Chemistry::addreax()
{
    std::string newname, style;
    ss >> newname;
    ss >> style;
    
    if (strcmp(style.c_str(), "simple") == 0){
        
        // record reaction name and typical arguments (gx calculator, DG caclulator, eq constant)
        Nreax++;
        rxnames.push_back(newname);
        std::string gxn,DGn;
        double Keqi, kii;
        ss >> gxn >> DGn >> Keqi >> kii;
        for (int i = 0; i<Ngx; i++) {
            if (strcmp(gxn.c_str(), gxnames[i].c_str())==0) rx_gxID.push_back(i);
        }
        for (int i = 0; i<NDG; i++) {
            if (strcmp(DGn.c_str(), DGnames[i].c_str())==0) rx_DGID.push_back(i);
        }
        if (rx_gxID.size()!=Nreax) {
            std::string msg = "ERROR: \""+gxn+"\" specificed for reaction \""+newname+"\" does not match any existing gammax calculator \n";
            error->errsimple(msg);
        }
        if (rx_DGID.size()!=Nreax) {
            std::string msg = "ERROR: \""+DGn+"\" specificed for reaction \""+newname+"\" does not match any existing DGx calculator \n";
            error->errsimple(msg);
        }
        Keq.push_back(Keqi);
        if (kii >= 0. && kii <= 1.) ki.push_back(kii);
        else {
            std::string msg = "ERROR: ki parameter specificed for reaction \""+newname+"\" must be between 0 and 1 \n";
            error->errsimple(msg);
        }
        
        
        // record list of stoichio changes associated with reaction, in both background and foreground
        bool bkg_on = false;
        bool fgd_on = false;
        bool solv_on = false;
        std::string entry, tmol,solvent;
        int nmol;
        std::vector<double> bvec,fvec;
        std::vector<int> ibvec,ifvec;
        std::vector<bool> isnotsolvec;
        bool comment_found = false;

        ss >> entry;  // IMPORTANT: you must read the entry before the operations so that while(ss) can recognize the end_of_stream value and stop reading at the right time
        while (ss && !comment_found){
            if (strcmp(entry.c_str(), "bkg") == 0) {
                bkg_on = true;
                fgd_on = false;
                solv_on = false;
                ss >> tmol;
            }
            else if (strcmp(entry.c_str(), "fgd") == 0) {
                bkg_on = false;
                fgd_on = true;
                solv_on=false;
                ss >> tmol;
            }
            else if (strcmp(entry.c_str(), "solvent") == 0) {
                bkg_on = false;
                fgd_on = false;
                solv_on=true;
                ss >> solvent;
            }
            else if (strncmp(entry.c_str(), "#", 1) == 0) {
                comment_found = true;
                bkg_on = false;
                fgd_on = false;
                solv_on = false;
            }
            else {
                tmol = entry;
            }
    
            
            if (bkg_on){
                bool bkgfound = false;
                for (int i=0; i<Nmol; i++) {
                    if (strcmp(tmol.c_str(), molnames[i].c_str())==0) {
                        ibvec.push_back(i);
                        bkgfound = true;
                        isnotsolvec.push_back(true);   //NB: if true, the molecule is not a solvent
                    }
                }
                if (!bkgfound) {
                    std::string msg = "ERROR: Molecule \""+tmol+"\" specificed for reaction \""+newname+"\" does not match any existing molecule \n";
                    error->errsimple(msg);
                }
                ss >> nmol;
                bvec.push_back(nmol);
            }
            
            if (fgd_on){
                bool fgdfound = false;
                for (int i=0; i<Nmol; i++) {
                    if (strcmp(tmol.c_str(), molnames[i].c_str())==0) {
                        ifvec.push_back(i);
                        fgdfound = true;
                    }
                }
                if (!fgdfound) {
                    std::string msg = "ERROR: Molecule \""+tmol+"\" specificed for reaction \""+newname+"\" does not match any existing molecule \n";
                    error->errsimple(msg);
                }
                ss >> nmol;
                fvec.push_back(nmol);
            }
            
            if (solv_on){
                // check in background vector created so far if solvent molecule type exists: if it does, turn its solvent flag on
                bool solvfound = false;
                for (int i=0; i<bvec.size(); i++) {
                    if (strcmp(solvent.c_str(), molnames[ibvec[i]].c_str())==0) {
                        isnotsolvec[i]=false;   //NB: if false , the molecule is a solvent so will be excluded from supersaturation
                        solvfound = true;
                    }
                }
                if (!solvfound) {
                    std::string msg = "ERROR: Molecule \""+solvent+"\" specificed as solvent for reaction \""+newname+"\" does not match any existing molecule in reaction's bkg defined so far\n";
                    error->errsimple(msg);
                }
            }
            
            
            ss >> entry; // IMPORTANT to keep this here and not at beginning of while(ss), otherwise end_of_stream is not processed at the right time
        }
        
        bkg_molID.push_back(ibvec);
        fgd_molID.push_back(ifvec);
        bkg_nmol.push_back(bvec);
        fgd_nmol.push_back(fvec);
        bkg_isnotsolv.push_back(isnotsolvec);
        
        
        // geometric changes to background and foreground induced by the reaction
        rx_dV_bkg.push_back(0.);
        for (int i=0; i<bvec.size(); i++)   rx_dV_bkg.back() += mol_vapp[ibvec[i]]*bvec[i];
        
        rx_dV_fgd.push_back(0.);
        for (int i=0; i<fvec.size(); i++)   rx_dV_fgd.back() += (mol_arad[ifvec[i]]*mol_acir[ifvec[i]]*mol_acir[ifvec[i]]) * fvec[i];

        rx_ar_min.push_back(0.);
        rx_ar_max.push_back(0.);
        rx_ar_cum.push_back(0.);
        double min = mol_arad[ifvec[0]];
        double max = mol_arad[ifvec[0]];
        double sum = 0.;
        for (int i=1; i<fvec.size(); i++){
            if (mol_arad[ifvec[i]] < min) min = mol_arad[ifvec[i]] ;
            if (mol_arad[ifvec[i]] > max) max = mol_arad[ifvec[i]] ;
            sum += mol_arad[ifvec[i]];
        }
        rx_ar_min.back() = min;
        rx_ar_max.back() = max;
        rx_ar_cum.back() = sum;
        rx_ar_avg.push_back(sum/((double)fvec.size()));
        rx_ar_avv.push_back(pow(rx_dV_fgd.back(),1./3.));
        
        rx_ac_min.push_back(0.);
        rx_ac_max.push_back(0.);
        rx_ac_cum.push_back(0.);
        min = mol_acir[ifvec[0]];
        max = mol_acir[ifvec[0]];
        sum = 0.;
        for (int i=1; i<fvec.size(); i++){
            if (mol_acir[ifvec[i]] < min) min = mol_acir[ifvec[i]] ;
            if (mol_acir[ifvec[i]] > max) max = mol_acir[ifvec[i]] ;
            sum += mol_acir[ifvec[i]];
        }
        rx_ac_min.back() = min;
        rx_ac_max.back() = max;
        rx_ac_cum.back() = sum;
        rx_ac_avg.push_back(sum/((double)fvec.size()));
        rx_ac_avv.push_back(pow(rx_dV_fgd.back(),1./3.));

        
        
        // print read simple reaction to screen to debug
        if (me==MASTER) {
            fprintf(screen," New simple reaction added: %s \n", rxnames[Nreax-1].c_str());
            fprintf(screen," Gammax calculator %s with ID %d \n", gxnames[rx_gxID[Nreax-1]].c_str() , rx_gxID[Nreax-1]);
            fprintf(screen," DGx calculator %s with ID %d \n", DGnames[rx_gxID[Nreax-1]].c_str(),rx_DGID[Nreax-1]);
            fprintf(screen," Equilibrium constant %f \n",Keq[Nreax-1]);
            fprintf(screen," Ki parameter %f \n",ki[Nreax-1]);
            fprintf(screen,"\n BKG: ");
            for (int i=0; i<(bkg_molID.back()).size(); i++) {
                fprintf(screen," %s %f",molnames[(bkg_molID.back())[i]].c_str(),(bkg_nmol.back())[i]);
            }
            fprintf(screen,"\n FGD: ");
            for (int i=0; i<(fgd_molID.back()).size(); i++) {
                fprintf(screen," %s %f",molnames[(fgd_molID.back())[i]].c_str(),(fgd_nmol.back())[i]);
            }
            fprintf(screen,"\n\n");
        }
        
        
    }
    else if (strcmp(style.c_str(), "chain") == 0){
        
        // add chain
        Nchain++;
        
        ch_dV_fgd.push_back(0.);
        ch_dV_bkg.push_back(0.);
        chnames.push_back(newname);
        ch_arac.push_back(1.);

        ss >> style;
        if (strcmp(style.c_str(), "parall") != 0  && strcmp(style.c_str(), "series") != 0) {
            std::string msg = "ERROR: Unknown style of chain: \""+newname+"\". Possible styles are parall or series, ibut \""+style+"\" was foud instead \n";
            error->errsimple(msg);
        }
        chstyle.push_back(style);
       
        std::vector<double> rdVvec_fgd,rdVvec_bkg;  // will be used to compute relative volume change due to each reaction in chain


        std::string entry, trx;
        int tnrx;
        bool comment_found = false;
        std::vector<int> irxvec,nrxvec;
        ss >> entry;  // IMPORTANT: you must read the entry before the operations so that while(ss) can recognize the end_of_stream value and stop reading at the right time
        while (ss && !comment_found){
            if (strncmp(entry.c_str(), "#", 1) == 0) comment_found = true;
            else if (strcmp(entry.c_str(), "arac") == 0) {
                double arac;
                ss >> arac;
                ch_arac[Nchain-1] = arac;
                ss >> entry;
            }
            else {
                trx = entry;
                ss >> tnrx;
                bool rxfound = false;
                for (int i=0; i<Nreax; i++) {
                    if (strcmp(trx.c_str(), rxnames[i].c_str()) == 0) {
                        irxvec.push_back(i);
                        ch_dV_fgd[Nchain-1] += ( rx_dV_fgd[i] * tnrx );
                        ch_dV_bkg[Nchain-1] += ( rx_dV_bkg[i] * tnrx );
                        
                        rdVvec_fgd.push_back(rx_dV_fgd[i] * tnrx);
                        rdVvec_bkg.push_back(rx_dV_bkg[i] * tnrx);
                        
                        rxfound = true;
                    }
                }
                if (!rxfound) {
                    std::string msg = "ERROR: Reaction \""+trx+"\" specificed for chain \""+newname+"\" does not match any existing reaction \n";
                    error->errsimple(msg);
                }
                nrxvec.push_back(tnrx);
                ss >> entry;
            }
        }
        ch_rxID.push_back(irxvec);
        ch_nrx.push_back(nrxvec);
        
        // compute relative volume change due to each reaction in chain
        for (int i=0; i<irxvec.size(); i++) {
            rdVvec_fgd[i] /= ch_dV_fgd[Nchain-1];
            rdVvec_bkg[i] /= ch_dV_bkg[Nchain-1];
        }
        ch_rdV_fgd.push_back(rdVvec_fgd);
        ch_rdV_bkg.push_back(rdVvec_bkg);
        
        
        // print to screen to debug
        if (me==MASTER) {
            fprintf(screen," New chain of reactions added: %s \n", chnames[Nchain-1].c_str());
            fprintf(screen," Reactions in %s are: ",chstyle[Nchain-1].c_str() );
            for (int i=0; i<(ch_rxID.back()).size(); i++) {
                fprintf(screen," %s %d",rxnames[(ch_rxID.back())[i]].c_str(),(ch_nrx.back())[i]);
            }
            fprintf(screen,"\n Relative fgd volume changes of reactions are: ");
            for (int i=0; i<(ch_rdV_fgd.back()).size(); i++) {
                fprintf(screen," %f",ch_rdV_fgd[Nchain-1][i]);
            }
            fprintf(screen,"\n\n");
        }
        
     
    }
    else {
        std::string msg = "ERROR: Unknown style of reaction or chain: \""+style+"\" \n";
        error->errsimple(msg);
    }
}
    


// -------------------------------------
// Add surface energy calculator
void Chemistry::addsen()
{
    Nsen++;
    
    std::string newname;
    ss >> newname;
    sennames.push_back(newname);
    
    ss >> newname;
    senstyle.push_back(newname);
    
    std::vector<double> cvec;
    if (strcmp(newname.c_str(), "const") == 0){
        double entry;
        ss >> entry;
        cvec.push_back(entry);
    }
    else if (strcmp(newname.c_str(), "tolman") == 0) {
        double e1,e2;
        ss >> e1 >> e2;
        cvec.push_back(e1);
        cvec.push_back(e2);
    }
    else if (strcmp(newname.c_str(), "lin") == 0) {
        double foo=1.;  // to be implemented..
    }
    else {
        std::string msg = "ERROR: Unknown style of gammax calculator: \""+newname+"\" \n";
        error->errsimple(msg);
    }
    
    sencoef.push_back(cvec);
    
    
    //just a check
    if (me==MASTER) fprintf(screen," New surface calculator added: %s Style = %s\n",sennames[Nsen-1].c_str(),senstyle[Nsen-1].c_str());
    if (me==MASTER) fprintf(screen," Current surface coeff matrix is: \n");
    for (int i=0; i<Nsen; i++) {
        for (int j=0; j<sencoef[i].size(); j++) {
            if (me==MASTER) fprintf(screen," %f ",sencoef[i][j]);
        }
        if (me==MASTER) fprintf(screen,"\n");
    }
}




// -------------------------------------
// Add mechanism
void Chemistry::addmech()
{
    Nmech++;
    
    std::string newname;
    ss >> newname;
    mechnames.push_back(newname);
    
    // read syle
    ss >> newname;
    if (strcmp(newname.c_str(), "allpar") != 0 && strcmp(newname.c_str(), "allser") != 0 && strcmp(newname.c_str(), "growth") != 0 && strcmp(newname.c_str(), "cnt") != 0) {
        std::string msg = "ERROR: Unknown style of mechanism \""+mechnames.back()+"\". Known styles are allpar, allser, growth, or cnt, instead \""+newname+"\" was found \n";
        error->errsimple(msg);
    }
    mechstyle.push_back(newname);
    
    // read mode
    ss >> newname;
    if (strcmp(newname.c_str(), "net") != 0 && strcmp(newname.c_str(), "straight") != 0) {
        std::string msg = "ERROR: Unknown mode of mechanism \""+mechnames.back()+"\". Known modes are net or straight, instead \""+newname+"\" was found \n";
        error->errsimple(msg);
    }
    mechmode.push_back(newname);
    
    
    // read name of reaction or chain
    ss >> newname;
    bool rxfound = false;
    for (int i = 0; i<Nreax; i++) {
        if (strcmp(newname.c_str(), rxnames[i].c_str()) == 0) {
            mechchain.push_back(false);
            mechrcID.push_back(i);
            rxfound = true;
        }
    }
    if (!rxfound){
        for (int i = 0; i<Nchain; i++) {
            if (strcmp(newname.c_str(), chnames[i].c_str()) == 0) {
                mechchain.push_back(true);
                mechrcID.push_back(i);
                rxfound = true;
            }
        }
    }
    if (!rxfound){
        std::string msg = "ERROR: Reaction or chain name \""+newname+"\" specified for mechanism\""+mechnames.back()+"\" does not correspond to any existing reactin or chain name \n";
        error->errsimple(msg);
    }

    
    // read surface energy calculator
    ss >> newname;
    bool senfound = false;
    for (int i = 0; i<Nsen; i++) {
        if (strcmp(newname.c_str(), sennames[i].c_str()) == 0) {
            mechsenID.push_back(i);
            senfound = true;
        }
    }
    if (!senfound){
        std::string msg = "ERROR: Reaction or chain name \""+newname+"\" specified for mechanism\""+mechnames.back()+"\" does not correspond to any existing reactin or chain name \n";
        error->errsimple(msg);
    }
    
    
    
    // read interaction energy scaling rule with particle size
    ss >> newname;
    if (strcmp(newname.c_str(), "int_no") != 0 && strcmp(newname.c_str(), "int_lin") != 0 && strcmp(newname.c_str(), "int_1lin") != 0 && strcmp(newname.c_str(), "int_2lin") != 0 && strcmp(newname.c_str(), "int_size") != 0 ) {
        std::string msg = "ERROR: Unknown energy scaling rule of  \""+mechnames.back()+"\". Known styles are int_no, int_lin, or int_size, or int_1Dlin, instead \""+newname+"\" was found \n";
        error->errsimple(msg);
    }
    mechinter.push_back(newname);
    
    
    //just a check
    if (me==MASTER) fprintf(screen," New mechanism added: %s Style = %s  Mode = %s\n",mechnames[Nmech-1].c_str(),mechstyle[Nmech-1].c_str(),mechmode[Nmech-1].c_str());
    if (me==MASTER) {
        if (mechchain.back())  fprintf(screen," Reaction is a CHAIN named : %s\n",chnames[mechrcID.back()].c_str());
        else fprintf(screen," Reaction is a SINGLE one named : %s\n",rxnames[mechrcID.back()].c_str());
    }
    if (me==MASTER) fprintf(screen," Interaction scaling: %s \n",mechinter[Nmech-1].c_str());
    if (me==MASTER) fprintf(screen,"\n");
}
    



// ---------------------------------------------------------------
// Compute DGx following the rules given by the user for the reaction taking place. This function is called when computing rates
double Chemistry::compDGx(int rxid)
{
    // for that reaction, find GDx calculator
    // execute rules of calculator
    double DGx = 0.;
    int DGid = rx_DGID[rxid];
    
    if (strcmp(DGstyle[DGid].c_str(),"const") == 0) DGx = DGcoef[DGid][0];
    else if (strcmp(DGstyle[DGid].c_str(),"test3coef") == 0) {
        //do nothing.. just a placeholder, but in this way you do not get an error
    }
    else {
        fprintf(screen,"\n ****************************************************** \n ****************************************************** chemistry.cpp -- ERROR: a reaction is trying to invoke a DG calculator whose sytle is not defined.. this should never happen, as it should be spotted already by the chemDB reader... something fishy is happening ******************************************************\n ******************************************************\n ******************************************************\n");
    }

    return DGx;  // the computed value of DGx is returned to the rate
}



// ---------------------------------------------------------------
// Compute gammax following the rules given by the user for the reaction taking place. This function is called when computing rates
double Chemistry::compgammax(int rxid)
{
    double gammax = 1.;
    int gID = rx_gxID[rxid];
    
    if (strcmp(gxstyle[gID].c_str(),"const") == 0) gammax = gxcoef[gID][0];
    else if (strcmp(gxstyle[gID].c_str(),"test3coef") == 0) {
        //do nothing.. just a placeholder, but in this way you do not get an error
    }
    else {
        fprintf(screen,"\n ****************************************************** \n ****************************************************** chemistry.cpp -- ERROR: a reaction is trying to invoke a Gamma* calculator whose sytle is not defined.. this should never happen, as it should be spotted already by the chemDB reader... something fishy is happening ******************************************************\n ******************************************************\n ******************************************************\n");
    }

    return gammax;  // the computed value of gammax is returned to the rate
    
}





// ---------------------------------------------------------------
// Compute surface energy following the rules given by the user for the reaction taking place. This function is called when computing rates
double Chemistry::compSen(int mid)
{
    double Sen = 1.;
    int sID = mechsenID[mid];
    
    if (strcmp(senstyle[sID].c_str(),"const") == 0) Sen = sencoef[sID][0];
    else if (strcmp(senstyle[sID].c_str(),"lin") == 0) {
        //do nothing.. just a placeholder, but in this way you do not get an error
    }
    else if (strcmp(senstyle[sID].c_str(),"tolman") == 0) {
        //do nothing.. just a placeholder, but in this way you do not get an error
    }
    else {
        fprintf(screen,"\n ****************************************************** \n ****************************************************** chemistry.cpp -- ERROR: a reaction is trying to invoke a Surface Energy calculator whose sytle is not defined.. this should never happen, as it should be spotted already by the chemDB reader... something fishy is happening ******************************************************\n ******************************************************\n ******************************************************\n");
    }
    
    return Sen;  // the computed value of surface energy is returned to the rate
    
}


    

// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Chemistry::printall()
{
	fprintf(screen,"\n---------ALL ABOUT chemistry----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
