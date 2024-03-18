#include "chemistry.h"
#include "memory.h"
#include "error.h"


//#include "universe.h"
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
    e0 = nullptr;
}



// ---------------------------------------------------------------
// Class destructor
Chemistry::~Chemistry()
{
    //delete lmp;
    
    if (!(e0 == nullptr)){
        free(e0[0]);
        free(e0);
        e0 = nullptr;
    }
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
            else if (strcmp(word.c_str(), "molecule_modify") == 0)   mol_modify();
            else if (strcmp(word.c_str(), "reaction_modify") == 0)   reax_modify();
            else if (strcmp(word.c_str(), "chain_modify") == 0)   ch_modify();
            else if (strcmp(word.c_str(), "mech_modify") == 0)   mech_modify();
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
    
    // compute volume of (solid) molecule
    mol_vm.push_back(mol_arad[Nmol-1] * mol_acir[Nmol-1] * mol_acir[Nmol-1]);
    
    mol_rcr0.push_back(size);
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
    if (me==MASTER) fprintf(screen," New molecule added: %s %f %f %f %f %f %f\n",molnames[Nmol-1].c_str(),mol_arad[Nmol-1],mol_acir[Nmol-1],mol_vapp[Nmol-1],mol_ahyd[Nmol-1],mol_z[Nmol-1], mol_vm[Nmol-1]);
    
    // assigning default valus to additional properties
    mol_Em.push_back(1.);
    mol_Gm.push_back(1.);
    mol_Us.push_back(0.);
    mol_Pr.push_back(0.);
    
    //just a check
  if (me==MASTER) fprintf(screen," Defaul additional propertis for the molecule: %f %f %f %f \n",mol_Em[Nmol-1],mol_Gm[Nmol-1],mol_Us[Nmol-1],mol_Pr[Nmol-1]);
  

}


// ---------------------------------------------------------------
// Add parameters to existing molecule
void Chemistry::mol_modify()
{
    // read molecule name
    std::string newname;
    ss >> newname;
    
    // find molecule ID
    int molID = -1;
    for (int i=0; i<molnames.size(); i++){
        if (strcmp(newname.c_str(), molnames[i].c_str()) == 0) molID = i;
    }
    if (molID == -1){
        std::string msg = "\nERROR: Unknown molecule invoked by molecule_modify command: "+newname+"\n\n";
        error->errsimple(msg);
    }
    
    // read through the rest of entries and, if recognized, record value for corresponding property
    while (ss >> newname){
        if (strcmp(newname.c_str(), "Em") == 0)  {
            ss >> mol_Em[molID];
            if (me==MASTER) fprintf(screen," Added property for molecule %s: Em = %f\n",molnames[molID].c_str(),mol_Em[molID]);
        }
        else if (strcmp(newname.c_str(), "Gm") == 0)  {
            ss >> mol_Gm[molID];
            if (me==MASTER) fprintf(screen," Added property for molecule %s: Gm = %f\n",molnames[molID].c_str(),mol_Gm[molID]);
        }
        else if (strcmp(newname.c_str(), "Us") == 0)  {
            ss >> mol_Us[molID];
            if (me==MASTER) fprintf(screen," Added property for molecule %s: Us = %f\n",molnames[molID].c_str(),mol_Us[molID]);
        }
        else if (strcmp(newname.c_str(), "Pr") == 0)  {
            ss >> mol_Pr[molID];
            if (me==MASTER) fprintf(screen," Added property for molecule %s: Pr = %f\n",molnames[molID].c_str(),mol_Pr[molID]);
        }
        else {
            std::string msg = "\nERROR: Unknown molecule property invoked by molecule_modify command for molecule  "+molnames[molID]+": "+newname+"\n\n";
            error->errsimple(msg);
        }
    }

}



// ---------------------------------------------------------------
// Add parameters to existing reaction
void Chemistry::reax_modify()
{
    // read reaction name
    std::string newname;
    ss >> newname;
    
    // find reaction ID
    int rxID = -1;
    for (int i=0; i<rxnames.size(); i++){
        if (strcmp(newname.c_str(), rxnames[i].c_str()) == 0) rxID = i;
    }
    if (rxID == -1){
        std::string msg = "\nERROR: Unknown reaction invoked by reaction_modify command: "+newname+"\n\n";
        error->errsimple(msg);
    }
    
    // read through the rest of entries and, if recognized, record value for corresponding property
    while (ss >> newname){
        if (strcmp(newname.c_str(), "Fk") == 0)  {
            ss >> rx_Fk[rxID];
            if (me==MASTER) fprintf(screen," Added property to reaction %s: Fk = %f\n",rxnames[rxID].c_str(),rx_Fk[rxID]);
        }
        else if (strcmp(newname.c_str(), "Uk") == 0)  {
            ss >> rx_Uk[rxID];
            if (me==MASTER) fprintf(screen," Added property to reaction %s: Uk = %f\n",rxnames[rxID].c_str(),rx_Uk[rxID]);
        }
        else {
            std::string msg = "\nERROR: Unknown reaction property invoked by reaction_modify command for reaction  "+rxnames[rxID]+": "+newname+"\n\n";
            error->errsimple(msg);
        }
    }

}


// ---------------------------------------------------------------
// Add parameters to existing chain
void Chemistry::ch_modify()
{
    // read chain name
    std::string newname;
    ss >> newname;
    
    // find chain ID
    int chID = -1;
    for (int i=0; i<chnames.size(); i++){
        if (strcmp(newname.c_str(), chnames[i].c_str()) == 0) chID = i;
    }
    if (chID == -1){
        std::string msg = "\nERROR: Unknown chain invoked by chain_modify command: "+newname+"\n\n";
        error->errsimple(msg);
    }
    
    // read through the rest of entries and, if recognized, record value for corresponding property
    while (ss >> newname){
        if (strcmp(newname.c_str(), "Fk") == 0)  {
            ss >> ch_Fk[chID];
            if (me==MASTER) fprintf(screen," Added property to chain %s: Fk = %f\n",chnames[chID].c_str(),ch_Fk[chID]);
        }
        else {
            std::string msg = "\nERROR: Unknown chain property invoked by chain_modify command for chain  "+chnames[chID]+": "+newname+"\n\n";
            error->errsimple(msg);
        }
    }

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
    
    // after the coefficients to compute DGx, the code expects the standard state concentration and dimensionality of the activated complex associated to DG
    
    ss >> newname;
    if (strcmp(newname.c_str(), "cx") != 0){
        std::string msg = "ERROR: In Chemistry database, expected cx after DGx coefficients for "+ DGnames.back()+", instad \""+newname+"\" was found \n";
        error->errsimple(msg);
    }
    else{
        double entry;
        ss >> entry;
        cx.push_back(entry);
    }
    
    ss >> newname;
    if (strcmp(newname.c_str(), "dim") != 0){
        std::string msg = "ERROR: In Chemistry database, expected dim after DGx coefficients for "+ DGnames.back()+", instad \""+newname+"\" was found \n";
        error->errsimple(msg);
    }
    else{
        double entry;
        ss >> entry;
        dim.push_back(entry);
    }
    
    
    //just a check
    if (me==MASTER) fprintf(screen," New DGx calculator added: %s Style = %s\n",DGnames[NDG-1].c_str(),DGstyle[NDG-1].c_str());
    if (me==MASTER) fprintf(screen," Current DGx coeff matrix , cx, and dim are: \n");
    for (int i=0; i<NDG; i++) {
        for (int j=0; j<DGcoef[i].size(); j++) {
            if (me==MASTER) fprintf(screen," %f ",DGcoef[i][j]);
        }
        if (me==MASTER) fprintf(screen,", cx = %e , dim = %f",cx[i],dim[i]);
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
        ki.push_back(kii);
        
        
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
        for (int i=0; i<bvec.size(); i++)   rx_dV_bkg.back() += mol_vapp[ibvec[i]] * bvec[i];
        
        rx_dV_fgd.push_back(0.);
        rx_dVp_fgd.push_back(0.);
        rx_dVt_fgd.push_back(0.);
        for (int i=0; i<fvec.size(); i++)  {
            // dV will be used to compute dV of chain steps, which are used to distribute interaction energy between steps in chain during fix_delete and fix_nucleate. Hence it does not include phase porosity
            rx_dV_fgd.back() += (mol_arad[ifvec[i]]*mol_acir[ifvec[i]]*mol_acir[ifvec[i]]) * fvec[i];
            
            // dVt will be used to compute molecular cross-sectional area of rates per unit surface. Hence must include phase porosity
            rx_dVt_fgd.back() += (mol_arad[ifvec[i]]*mol_acir[ifvec[i]]*mol_acir[ifvec[i]]) * fvec[i] * mol_rcr0[ifvec[i]] * (1.+mol_Pr[ifvec[i]]);
            // dVp = sam as dV above, but including phase porosity (to be used to compute number of units to dissolve or nucleate in a particle (nrv in fix delete and nucleate)
            rx_dVp_fgd.back() += (mol_arad[ifvec[i]]*mol_acir[ifvec[i]]*mol_acir[ifvec[i]]) * fvec[i] * (1.+mol_Pr[ifvec[i]]);
        }

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

        rx_Fk.push_back(1.);     // Initialising kink fraction for new reaction. This can be modified later via the reaction_modify input command
        
        rx_Uk.push_back(0.);     // Default energy of reaction in kink position in zero
        
        // print read simple reaction to screen to debug
        if (me==MASTER) {
            fprintf(screen," New simple reaction added: %s \n", rxnames[Nreax-1].c_str());
            fprintf(screen," Gammax calculator %s with ID %d \n", gxnames[rx_gxID[Nreax-1]].c_str() , rx_gxID[Nreax-1]);
            fprintf(screen," DGx calculator %s with ID %d \n", DGnames[rx_gxID[Nreax-1]].c_str(),rx_DGID[Nreax-1]);
            fprintf(screen," Equilibrium constant %f \n",Keq[Nreax-1]);
            fprintf(screen,"\n BKG: ");
            for (int i=0; i<(bkg_molID.back()).size(); i++) {
                fprintf(screen," %s %f",molnames[(bkg_molID.back())[i]].c_str(),(bkg_nmol.back())[i]);
            }
            fprintf(screen,"\n FGD: ");
            for (int i=0; i<(fgd_molID.back()).size(); i++) {
                fprintf(screen," %s %f",molnames[(fgd_molID.back())[i]].c_str(),(fgd_nmol.back())[i]);
            }
            fprintf(screen," FGD volume change without (%e) and with (%e) pores\n", rx_dV_fgd[Nreax-1], rx_dVp_fgd[Nreax-1]);
            fprintf(screen,"\n\n");
        }
        
        
    }
    else if (strcmp(style.c_str(), "chain") == 0){
        
        // add chain
        Nchain++;
        
        ch_dV_fgd.push_back(0.);
        ch_dVp_fgd.push_back(0.);
        ch_dV_bkg.push_back(0.);
        chnames.push_back(newname);
        ch_arac.push_back(1.);

        ss >> style;
        if (strcmp(style.c_str(), "parall") != 0  && strcmp(style.c_str(), "series") != 0) {
            std::string msg = "ERROR: Unknown style of chain: \""+newname+"\". Possible styles are parall or series, ibut \""+style+"\" was foud instead \n";
            error->errsimple(msg);
        }
        chstyle.push_back(style);
       
        std::vector<double> rdVvec_fgd, rdVvec_bkg, rdVpvec_fgd;  // will be used to compute relative volume change due to each reaction in chain


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
                        ch_dVp_fgd[Nchain-1] += ( rx_dVp_fgd[i] * tnrx );
                        ch_dV_bkg[Nchain-1] += ( rx_dV_bkg[i] * tnrx );
                        
                        rdVvec_fgd.push_back(rx_dV_fgd[i] * tnrx);
                        rdVpvec_fgd.push_back(rx_dVp_fgd[i] * tnrx);
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
        
        ch_Fk.push_back(1.);    // Initialising kink fraction for the chain. This can be changed later using chain_modify in the chemDB input file
        
        // compute relative volume change due to each reaction in chain
        for (int i=0; i<irxvec.size(); i++) {
            rdVvec_fgd[i] /= ch_dV_fgd[Nchain-1];
            rdVpvec_fgd[i] /= ch_dVp_fgd[Nchain-1];
            rdVvec_bkg[i] /= ch_dV_bkg[Nchain-1];
        }
        ch_rdV_fgd.push_back(rdVvec_fgd);
        ch_rdVp_fgd.push_back(rdVpvec_fgd);
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
                fprintf(screen," %f, %f",ch_rdV_fgd[Nchain-1][i],ch_rdVp_fgd[Nchain-1][i]);
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
    if (strcmp(newname.c_str(), "allpar") != 0 && strcmp(newname.c_str(), "allser") != 0 && strcmp(newname.c_str(), "growth") != 0 && strcmp(newname.c_str(), "cnt") != 0 && strcmp(newname.c_str(), "micro") != 0) {
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
    if (strcmp(newname.c_str(), "int_no") != 0 && strcmp(newname.c_str(), "int_lin") != 0 && strcmp(newname.c_str(), "int_1lin") != 0 && strcmp(newname.c_str(), "int_2lin") != 0 && strcmp(newname.c_str(), "int_size") != 0 && strcmp(newname.c_str(), "int_uni") != 0) {
        std::string msg = "ERROR: Unknown energy scaling rule of  \""+mechnames.back()+"\". Known styles are int_no, int_lin, or int_size, or int_1Dlin, instead \""+newname+"\" was found \n";
        error->errsimple(msg);
    }
    mechinter.push_back(newname);
    
    
    // default rate style is TST (with excess enthalpy correction, chi parameter, etc as usual in MASKE)
    mechrate.push_back("TST");
    
    
    std::vector<std::string> tempvec;
    int rows = 0;   // number of raws in e0, ef and gij matrices. This equals the largest particle type among user-listed real and trial types in input file
    
    if(strcmp(mechstyle.back().c_str(),"micro")==0){
        std::string param;
        std::string msg = "ERROR: Mechanism \""+mechnames.back()+"\" is of style \"micro\" which requires coverage_calculator_style and associated filename at the end of the \"mech\" command in chemDB\n";
        
        if(ss >> param) tempvec.push_back(param);   // coverage calculator style name
        else  error->errsimple(msg);
        
        if (strcmp(param.c_str(),"pair") != 0){
            msg = "ERROR: Unknown coverage calculator style \""+param+"\". For now, only style \"pair\" has been implemented. Future styles may be \"bond\", \"ellipse\", etc.\n";
            error->errsimple(msg);
        }
        
        if(ss >> param) tempvec.push_back(param);   // filename for coverage calculator
        else error->errsimple(msg);
        
        // copy content of coverage calculator file into e0, ef and gij arrays here...
        
        // ... first we allocate the arrays
        for (int i =0; i<msk->Rtypes.size(); i++){
            if (msk->Rtypes[i] > rows) rows = msk->Rtypes[i];
        }
        for (int i =0; i<msk->Ttypes.size(); i++){
            if (msk->Ttypes[i] > rows) rows = msk->Ttypes[i];
        }
        
        if (rows == 0){
            msg = "ERROR: The matrices associated to file \""+param+"\" have size = 0. In your input script, make sure to list the real and trial types before reading the ChemDB file\n";
            error->errsimple(msg);
        }
        
        int nbytes = ((int) sizeof(double)) * rows * rows;
        double *data = (double *) malloc(nbytes);
        double *data1 = (double *) malloc(nbytes);
        double *data2 = (double *) malloc(nbytes);
        nbytes = ((int) sizeof(double *)) * rows;
        e0 = (double **) malloc(nbytes);
        ef = (double **) malloc(nbytes);
        gij = (double **) malloc(nbytes);
        
        int n = 0;
        for (int i = 0; i < rows ; i++) {
            e0[i] = &data[n];
            ef[i] = &data1[n];
            gij[i] = &data2[n];
            n += rows;
        }
        
        // ... then we pre-populate them with "-1"
        for (int i=0; i<rows; i++){
            for (int j=0; j<rows; j++){
                e0[i][j]=-1.;
                ef[i][j]=-1.;
                gij[i][j]=-1.;
            }
        }
        
        
        // ... finally we go thorugh the coverage file and import the correct values...
        std::ifstream inPfile(param.c_str());
        if (!inPfile.is_open()) {
            std::string msg = "ERROR: Additional parameters file for micro mechanism not found: \""+param+"\"";
            error->errsimple(msg);
        }
        else {
            fprintf(screen,"PROC %d: Reading additional parameters file %s\n",me,param.c_str());

            while (!inPfile.eof()) {
                std::stringstream iss;
                std::string line, word;
                std::getline (inPfile, line);
                iss.clear();
                iss.str(line);
                fprintf(screen,"PROC %d: Reading line: %s\n",me,line.c_str());
                int it = -1;
                int jt = -1;
                int rcol = 0;
                bool errflag = false;
                if (iss >> word){
                    rcol = 1;
                    if (strncmp(word.c_str(), "#", 1) == 0) int foo = 1; // Do nothing. Lines starting with # are comments
                    else {
                        // read i type
                        bool istar = false;
                        if (strcmp(word.c_str(), "*") == 0) istar = true;
                        else {
                            it = atoi(word.c_str());
                            if (it < 1 || it > rows){
                                std::string msg = "\nERROR: In the additional parameters file (\""+param+"\") you have entered, in column 1, a type that is either <= 0 or greater than the largest type listed among real and trial types in the input script\n\n";
                                error->errsimple(msg);
                            }
                        }
                        
                        // read j type
                        bool jstar = false;
                        if (iss >> word){
                            if (strcmp(word.c_str(), "*") == 0) jstar = true;
                            else {
                                jt = atoi(word.c_str());
                                if (jt < 1 || jt > rows){
                                    std::string msg = "\nERROR: In the additional parameters file (\""+param+"\") you have entered, in column 1, a type that is either <= 0 or greater than the largest type listed among real and trial types in the input script\n\n";
                                    error->errsimple(msg);
                                }
                            }
                        }
                        else errflag = true;
                        
                        
                        
                        // record e0 value(s)
                        // NB: "-1" everywhere because LAMMPS types start from 1 and arrays in C++ start from 0
                        double val;
                        if (iss >> val){
                            if (!istar && !jstar){
                                e0[it-1][jt-1] = val;
                                e0[jt-1][it-1] = val;
                            }
                            else if(istar && !jstar){
                                for (int i=1; i<rows+1; i++){
                                    e0[i-1][jt-1] = val;
                                    e0[jt-1][i-1] = val;
                                }
                            }
                            else if(jstar && !istar){
                                for (int i=1; i<rows+1; i++){
                                    e0[i-1][it-1] = val;
                                    e0[it-1][i-1] = val;
                                }
                            }
                            else if(jstar && istar){
                                for (int i=1; i<rows+1; i++){
                                    for (int j=1; j<rows+1; j++){
                                        e0[i-1][j-1] = val;
                                    }
                                }
                            }
                        }
                        else errflag = true;
                        
                        
                        // record ef value(s)
                        // NB: "-1" everywhere because LAMMPS types start from 1 and arrays in C++ start from 0
                        if (iss >> val){
                            if (!istar && !jstar){
                                ef[it-1][jt-1] = val;
                                ef[jt-1][it-1] = val;
                            }
                            else if(istar && !jstar){
                                for (int i=1; i<rows+1; i++){
                                    ef[i-1][jt-1] = val;
                                    ef[jt-1][i-1] = val;
                                }
                            }
                            else if(jstar && !istar){
                                for (int i=1; i<rows+1; i++){
                                    ef[i-1][it-1] = val;
                                    ef[it-1][i-1] = val;
                                }
                            }
                            else if(jstar && istar){
                                for (int i=1; i<rows+1; i++){
                                    for (int j=1; j<rows+1; j++){
                                        ef[i-1][j-1] = val;
                                    }
                                }
                            }
                        }
                        else errflag = true;
                        
                        
                        // record gij value(s)
                        // NB: "-1" everywhere because LAMMPS types start from 1 and arrays in C++ start from 0
                        if (iss >> val){
                            if (!istar && !jstar){
                                gij[it-1][jt-1] = val;
                                gij[jt-1][it-1] = val;
                            }
                            else if(istar && !jstar){
                                for (int i=1; i<rows+1; i++){
                                    gij[i-1][jt-1] = val;
                                    gij[jt-1][i-1] = val;
                                }
                            }
                            else if(jstar && !istar){
                                for (int i=1; i<rows+1; i++){
                                    gij[i-1][it-1] = val;
                                    gij[it-1][i-1] = val;
                                }
                            }
                            else if(jstar && istar){
                                for (int i=1; i<rows+1; i++){
                                    for (int j=1; j<rows+1; j++){
                                        gij[i-1][j-1] = val;
                                    }
                                }
                            }
                        }
                        else errflag = true;


                        
                        if (errflag){
                            std::string msg = "\nERROR: Each row in the additional parameters file (\""+param+"\") must have 5 entries. At least one row has fewer - double check the file\n\n";
                            error->errsimple(msg);
                        }
                    }
                }
                    
                    
                
            }
        }
        
        
        //fprintf(screen,"DEBUG 0: PROC %d param value is %s",me,param.c_str());
        msg = "ERROR: Mechanism style requires an interaction energy calculator style \n";
        if(ss >> param) tempvec.push_back(param);   // interaction energy calculator style ("energy" or "stress")
        else error->errsimple(msg);
        //fprintf(screen,"DEBUG 1: PROC %d param value is %s",me,param.c_str());
        //sleep (1);
    
        if (strcmp(param.c_str(),"energy") != 0){
            msg = "ERROR: Unknown energy calculator style \""+param+"\". For now, only style \"energy\" has been implemented. Future style may be \"stress\".\n";
            error->errsimple(msg);
        }

        if(ss >> param) tempvec.push_back(param);   // limiting coverage fraction
        else{
            msg = "ERROR: Mechanism style requires a limiting coverage fraction \n";
            error->errsimple(msg);
        }
    
        

    }
    mechpar.push_back(tempvec);
    
    
    //just a check
    if (me==MASTER) fprintf(screen," New mechanism added: %s Style = %s  Mode = %s\n",mechnames[Nmech-1].c_str(),mechstyle[Nmech-1].c_str(),mechmode[Nmech-1].c_str());
    if (me==MASTER) {
        if (mechchain.back())  fprintf(screen," Reaction is a CHAIN named : %s\n",chnames[mechrcID.back()].c_str());
        else fprintf(screen," Reaction is a SINGLE one named : %s\n",rxnames[mechrcID.back()].c_str());
    }
    if (me==MASTER) fprintf(screen," Interaction scaling: %s \n",mechinter[Nmech-1].c_str());
    if (me==MASTER) fprintf(screen," Rate Style: %s \n",mechrate[Nmech-1].c_str());
    if (me==MASTER) fprintf(screen," Additional parameters:");
    if (me==MASTER) {
        for (int pp=0;pp<mechpar[Nmech-1].size();pp++){
            fprintf(screen," %s",mechpar[Nmech-1][pp].c_str());
        }
        fprintf(screen,"\n\n");
    }
    if (me==MASTER){
        if(strcmp(mechstyle.back().c_str(),"micro")==0){
            if (strcmp(mechpar[Nmech-1][0].c_str(),"pair") == 0){
                fprintf(screen," Values in e0 matrix from file %s \n",(mechpar[Nmech-1].back()).c_str());
                for (int i=0; i<rows; i++){
                    for (int j=0; j<rows; j++){
                        fprintf(screen,"%e ",e0[i][j]);
                    }
                    fprintf(screen,"\n");
                }
                
                fprintf(screen," Values in ef matrix from file %s\n",(mechpar[Nmech-1].back()).c_str());
                for (int i=0;i<rows;i++){
                    for (int j=0; j<rows; j++){
                        fprintf(screen,"%e ",ef[i][j]);
                    }
                    fprintf(screen,"\n");
                }
                
                fprintf(screen," Values in gij matrix from file %s\n",(mechpar[Nmech-1].back()).c_str());
                for (int i=0;i<rows;i++){
                    for (int j=0; j<rows; j++){
                        fprintf(screen,"%e ",gij[i][j]);
                    }
                    fprintf(screen,"\n");
                }
            }
        }
    }
    if (me==MASTER) fprintf(screen,"\n\n");
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
// Modify parameters of existing mechanism
void Chemistry::mech_modify()
{
    // read mechanism name
    std::string newname;
    ss >> newname;
    
    // find mecvhanism ID
    int mechID = -1;
    for (int i=0; i<mechnames.size(); i++){
        if (strcmp(newname.c_str(), mechnames[i].c_str()) == 0) mechID = i;
    }
    if (mechID == -1){
        std::string msg = "\nERROR: Unknown mechanism invoked by molecule_modify command: "+newname+"\n\n";
        error->errsimple(msg);
    }
    
    // read through the rest of entries and, if recognized, record value for corresponding property
    while (ss >> newname){
        if (strcmp(newname.c_str(), "rate") == 0)  {
            ss >> mechrate[mechID];
            if (me==MASTER) fprintf(screen," Modified rate style for mechanism: %s %s\n",mechnames[mechID].c_str(),mechrate[mechID].c_str());
        }
        else {
            std::string msg = "\nERROR: Unknown mechanism property invoked by mech_modify command in ChemDB, for mechanism  "+mechnames[mechID]+": "+newname+"\n\n";
            error->errsimple(msg);
        }
    }

}



    

// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Chemistry::printall()
{
	fprintf(screen,"\n---------ALL ABOUT chemistry----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
