#include "solution.h"
#include "chemistry.h"
#include "error.h"
#include "output.h"
#include "lammpsIO.h"
#include "universe.h"
#include "fix.h"

/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

*/

#include <string.h>

using namespace MASKE_NS;


Solution::Solution(MASKE *maske) : Pointers(maske)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    soltype = "uniform";
    dV = 0.;
    voidV = 0.;
    dVvoidV = 0.;
    unitC = -1.;
}



// ---------------------------------------------------------------
// Class destructor
Solution::~Solution()
{
    //delete lmp;
}



// ---------------------------------------------------------------
// add a molecule from input file
void Solution::addmol(std::string molname , double conci)
{
    molins.push_back(molname);
    
    bool found_mol=false;
    for (int i=0; i < chem->Nmol; i++) {
        if ( strcmp(molname.c_str(),(chem->molnames[i]).c_str())==0 ) {
            chem->mol_cins[i] = conci;
            found_mol = true;
        }
    }
    
    if (!found_mol) {
        std::string msg = "ERROR: Unknown molecule invoked by the sol_start command: \""+molname+"\" \n";
        error->errsimple(msg);
    }
    
    // recording concentrations from input file
    conc.push_back(conci);
    concdV.push_back(conci);
}



// ---------------------------------------------------------------
// function to compute the number of molecules in solution in box and in dV
void Solution::computeNmol(void)
{
    // processors with active lammps read box volume and particle radii from lammps and compute packing fraction
    if (lammpsIO->lammps_active){
        
        std::string tolmp;
        //tolmp = "compute tempRAD all property/atom radius";
        //lammpsIO->lammpsdo(tolmp);
        tolmp = "variable Bvol equal vol";
        lammpsIO->lammpsdo(tolmp);
        tolmp = "compute rad all property/atom radius"; // change all to exclude bacteria
        lammpsIO->lammpsdo(tolmp);
        tolmp = "variable sva atom 4./3.*PI*c_rad*c_rad*c_rad";
        lammpsIO->lammpsdo(tolmp);
        tolmp = "compute solidVc all reduce sum v_sva"; // change all to exclude bacteria
        lammpsIO->lammpsdo(tolmp);
        tolmp = "variable solidV equal c_solidVc";
        lammpsIO->lammpsdo(tolmp);
        tolmp = "variable packf equal v_solidV/vol";
        lammpsIO->lammpsdo(tolmp);

        tolmp = lammpsIO->lmpThSt + " c_solidVc v_solidV v_packf";
        lammpsIO->lammpsdo(tolmp);
        lammpsIO->lammpsdo("run 0");
        lammpsIO->lammpsdo(tolmp);
        
       
        int natoms  = static_cast<int> (lammpsIO->lmp->atom->natoms);
        BoxV = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"Bvol",0));

        //double *aR;
        //aR = new double[natoms];
        //aR = ((double *) lammps_extract_compute(lammpsIO->lmp,(char *) "tempRAD",1,1));
        
        //SolidV = 0.;
        SolidV = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)"solidV",0));
        
        /*for (int i =0; i<natoms; i++) {
            fprintf(screen,"\n Proc %d ..... R[%d] = %f \n",me,i,aR[i]);
            SolidV += 4. / 3. * M_PI * aR[i] * aR[i] * aR[i];   // ATTENTION: ONLY TRUE FOR SPHERICAL PARTICLES!
        }*/
        
       /* sleep(2);
        fprintf(screen,"\n Proc %d is OK",me);
        MPI_Barrier(MPI_COMM_WORLD);
        */
        PackF = SolidV/BoxV;
        
        //tolmp = "uncompute tempRAD";
        //lammpsIO->lammpsdo(tolmp);
        tolmp = "variable Bvol delete"; lammpsIO->lammpsdo(tolmp);
        tolmp = "uncompute rad";        lammpsIO->lammpsdo(tolmp);
        tolmp = "variable sva delete"; lammpsIO->lammpsdo(tolmp);
        tolmp = "uncompute solidVc";        lammpsIO->lammpsdo(tolmp);
        tolmp = "variable solidV delete"; lammpsIO->lammpsdo(tolmp);
        tolmp = "variable packf delete"; lammpsIO->lammpsdo(tolmp);
        lammpsIO->lammpsdo(lammpsIO->lmpThSt);    //restore original thermo
    }
    
  
    
    
    
    // sumbaster of first-lammps-subcomm sends packF and BoxV to all other processors
    if (me == universe->flampID) {
        for (int dest=0; dest<universe->nprocs; dest++) {
            if (dest!=me) {
                MPI_Send(&BoxV, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
                MPI_Send(&PackF, 1, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
            }
        }
    }
    else {
        int source = universe->flampID;
        MPI_Recv(&BoxV, 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&PackF, 1, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
    }
    // try
    // all procs compute volume of solution in box and in dV, and number if molecules depending on concentrations
    SVol = BoxV * (1.- PackF) - voidV;
    dVSVol = dV - dVvoidV;
    
    // number of molecules of each species in solution in box and in dV  (for now I assume same starting concentration in box and dV... to be changed in future)
    for (int i =0; i<conc.size(); i++){
        nmol.push_back(SVol * conc[i] * nAvo * unitC);
        for (int j=0; j < chem->Nmol; j++) {
            if ( strcmp(molins[i].c_str(),(chem->molnames[j]).c_str())==0 ) {
                chem->mol_nins[j] = nmol.back();
                chem->mol_cindV[j] = conc[i];
                chem->mol_nindV[j] = nmol.back()/SVol*dVSVol;
            }
        }
    }
    
    if (me == MASTER) {
        fprintf(screen, "\nAvogadro number is %e\nSolution volume is %f\n Packing fraction is %f\n Unit conversion vol is %e\n dV is %f\n dVvoidV is %f\n voidV is %f\n DH_A is %f\n",nAvo,SVol,PackF,unitC,dV,dVvoidV,voidV,DH_A);
    }
    
    
}





// ---------------------------------------------------------------
// compute activity product using the approach in PHREEQC (Truesdell-Jones equation if Debye Huckel parameters a and b in molecules in chemDB are >= 0. If they are < 0, uses Davies for charged molecules and Setschenow equation (log10 gamma = 0.1 Istrength) for uncharged ones
double Solution::compQ(int rxid, std::string type, std::string sol_in_style, std::string sol_in_UL)
{
    //double beta = 1.;
    double Q = 1.;   // activity product
    
    if (strcmp(sol_in_UL.c_str(),"uniform")==0   &&  strcmp(sol_in_style.c_str(),"fixed")==0  ) {
        // in future, I should be able to treat solutions that change during reaction steps in fixKMC. This will require update of virtual concentration vectors
        
        // a vector reading the concentrations of all molecules in solution  (PROBLEMATIC IF some charged molecules are considered explicitly as particles instead)
        std::vector<double> Vcins;
        for (int i=0; i<chem->Nmol; i++) {
            Vcins.push_back(chem->mol_cins[i] );
            if (Vcins[i]<0.) Vcins[i]=0.;
        }
        
        // compute ionic strength  (can be moved in update solution..)
        double Istr = 0.;
        for (int i=0; i<chem->Nmol; i++) {
            Istr +=  Vcins[i] * chem->mol_z[i] * chem->mol_z[i];
        }
        Istr /= 2.;
        double Isqrt = sqrt(Istr);
        
        std::string msg;
        msg = "";
        if (msk->wplog) {
            msg += "rxid ";
            std::ostringstream ss;    ss << rxid;   msg += ss.str(); ss.str("");   ss.clear();
            msg += ", Istr ";
            ss << Istr;   msg += ss.str();
        }
        
        
        
        if (strcmp(type.c_str(),"reac")==0) {
            //  go through the vector of bkg changes in chem for the current reaction
            for (int i=0; i<(chem->bkg_molID[rxid]).size(); i++) {
                // if the molecule involved in the reaction is not a solvent, and if it is consumed (thus negative number involved)
                if (chem->bkg_isnotsolv[rxid][i] && chem->bkg_nmol[rxid][i]<0) {
                    
                    double gam;
                    int molID = chem->bkg_molID[rxid][i];
                    
                    if (chem->mol_ahyd[molID] < 0. || chem->mol_bDH[molID] < 0.){
                        // if a Debye Huckel constants is input as negative, use Davies equation without them...
                        if (chem->mol_z[molID] != 0 ){
                            gam = pow( 10. , - DH_A * chem->mol_z[molID] * chem->mol_z[molID] * (Isqrt / ( 1.+ Isqrt) - 0.3 * Istr));
                        } else{
                            gam = pow( 10. , 0.1 * Istr);
                        }
                    }
                    else{
                        // ... otherwise use Debye Huckel equation with added explicit term on inonic strength for uncharged molecules and concentrated electrolites
                        gam = pow( 10. , - chem->mol_z[molID] * chem->mol_z[molID] * (DH_A * Isqrt) / ( 1.+ DH_B * chem->mol_ahyd[molID] * Isqrt) + chem->mol_bDH[molID]*Istr);
                    }
                    
                    // the below is   ( gamma * conc ) ^  stoichio
                    Q *= pow(gam * Vcins[molID] , (double) -chem->bkg_nmol[rxid][i]);
                    
                    if (msk->wplog) {
                        msg += ", mol ";
                        std::ostringstream ss;    ss << i;   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", molID ";
                        ss << molID;   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", virtual_mol_cins ";
                        ss << Vcins[molID];   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", mol_stoichio ";
                        ss << chem->bkg_nmol[rxid][i];   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", gam ";
                        ss << gam;   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", Q_reac ";
                        ss << Q;   msg += ss.str();
                        output->toplog(msg);
                    }
                }
            }
        }
        else if(strcmp(type.c_str(),"prod")==0){
            //  go through the vector of bkg changes in chem for the current reaction
            for (int i=0; i<(chem->bkg_molID[rxid]).size(); i++) {
                // if the molecule involved in the reaction is not a solvent, and if it is produced (thus positive number involved)
                if (chem->bkg_isnotsolv[rxid][i] && chem->bkg_nmol[rxid][i]>0) {
                    
                    double gam;
                    int molID = chem->bkg_molID[rxid][i];
                    
                     if (chem->mol_ahyd[molID] < 0. || chem->mol_bDH[molID] < 0.){
                         // if a Debye Huckel constants is input as negative, use Davies equation without them...
                         if (chem->mol_z[molID] != 0 ){
                             gam = pow( 10. , - DH_A * chem->mol_z[molID] * chem->mol_z[molID] * (Isqrt / ( 1.+ Isqrt) - 0.3 * Istr));
                        } else{
                            gam = pow( 10. , 0.1 * Istr);
                        }
                    }
                    else{
                        // ... otherwise use Debye Huckel equation with added explicit term on inonic strength for uncharged molecules and concentrated electrolites
                        gam = pow( 10. , - chem->mol_z[molID] * chem->mol_z[molID] * (DH_A * Isqrt) / ( 1.+ DH_B * chem->mol_ahyd[molID] * Isqrt) + chem->mol_bDH[molID]*Istr);
                    }
                    
                    // the below is   ( gamma * conc ) ^  stoichio
                    Q *= pow(gam * Vcins[molID] , (double) chem->bkg_nmol[rxid][i]);
                    
                    if (msk->wplog) {
                        msg += ", mol ";
                        std::ostringstream ss;    ss << i;   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", molID ";
                        ss << molID;   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", virtual_mol_cins ";
                        ss << Vcins[molID];   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", mol_stoichio ";
                        ss << chem->bkg_nmol[rxid][i];   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", gam ";
                        ss << gam;   msg += ss.str(); ss.str("");   ss.clear();
                        msg += ", Q_prod ";
                        ss << Q;   msg += ss.str();
                        output->toplog(msg);
                    }
                }
            }
        }
        
      
        
        if (msk->wplog) {
            if (strcmp(type.c_str(),"prod")==0)         msg += ", Qprod ";
            else if (strcmp(type.c_str(),"reac")==0)    msg += ", Qreac ";
            
            std::ostringstream ss;    ss << Q;   msg += ss.str();
        }
    }
    
    return Q;
}






// ---------------------------------------------------------------
// compute beta using the Truesdell or Davies and Setschenow equations
//   OBSOLETE as of 2021-02-05.... To be deleted when rationalising this class
double Solution::compbeta(int rxid, bool net, std::string sol_in_style, std::string sol_in_UL)
{
    double beta = 1.;
    double Q = 1.;   // activity product
    
    if (strcmp(sol_in_UL.c_str(),"uniform")==0   &&  strcmp(sol_in_style.c_str(),"fixed")==0  ) {
        // in future, I should be able to treat solutions that change during reaction steps in fixKMC. This will require update of virtual concentration vectors
        
        // a vector reading the concentrations of molecules in solution
        std::vector<double> Vcins;
        for (int i=0; i<chem->Nmol; i++) {
            Vcins.push_back(chem->mol_cins[i] );
            if (Vcins[i]<0.) Vcins[i]=0.;
        }
        
        // compute ionic strength
        double Istr = 0.;
        for (int i=0; i<chem->Nmol; i++) {
            Istr +=  Vcins[i] * chem->mol_z[i] * chem->mol_z[i];
        }
        Istr /= 2.;
        double Isqrt = sqrt(Istr);
        
        std::string msg;
        msg = "";
        if (msk->wplog) {
            msg += "rxid ";
            std::ostringstream ss;    ss << rxid;   msg += ss.str(); ss.str("");   ss.clear();
            msg += ", Istr ";
            ss << Istr;   msg += ss.str();
        }
        
        
        //  go through the vector of bkg changes in chem for the current reaction
        //   beta = prod (( gamma * conc)^stoichio_coeff)    / Keq
        for (int i=0; i<(chem->bkg_molID[rxid]).size(); i++) {

            // if the molecule involved in the reaction is not a solvent
            if (chem->bkg_isnotsolv[rxid][i]) {
                double gam;
                int molID = chem->bkg_molID[rxid][i];
                
                if (chem->mol_ahyd[molID] < 0. || chem->mol_bDH[molID] < 0.){
                     // if a Debye Huckel constants is input as negative, use Davies equation without them...
                     if (chem->mol_z[molID] != 0 ){
                         gam = pow( 10. , - DH_A * chem->mol_z[molID] * chem->mol_z[molID] * (Isqrt / ( 1.+ Isqrt) - 0.3 * Istr));
                    } else{
                        gam = pow( 10. , 0.1 * Istr);
                    }
                }
                else{
                    // ... otherwise use Debye Huckel equation with added explicit term on inonic strength for uncharged molecules and concentrated electrolites
                    gam = pow( 10. , - chem->mol_z[molID] * chem->mol_z[molID] * (DH_A * Isqrt) / ( 1.+ DH_B * chem->mol_ahyd[molID] * Isqrt) + chem->mol_bDH[molID]*Istr);
                }
                
                // the below is   ( gamma * conc ) ^  stoichio  
                Q *= pow(gam * Vcins[molID] , (double) chem->bkg_nmol[rxid][i]);
                
                if (msk->wplog) {
                    msg += ", mol ";
                    std::ostringstream ss;    ss << i;   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", molID ";
                    ss << molID;   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", virtual_mol_cins ";
                    ss << Vcins[molID];   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", mol_stoichio ";
                    ss << chem->bkg_nmol[rxid][i];   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", gam ";
                    ss << gam;   msg += ss.str(); ss.str("");   ss.clear();
                    msg += ", Q_temp ";
                    ss << Q;   msg += ss.str();
                }
            }
            
            
            
        }
        
        if (msk->wplog) {
            msg += ", Qprod ";
            std::ostringstream ss;    ss << Q;   msg += ss.str();
        }
        // if net = 1 reaction is backward, hence one should use Keq^-1 from user input
        if (net) beta = Q/chem->Keq[rxid];
        else beta = chem->Keq[rxid]/Q;
        
        if (msk->wplog) {
            msg += ", beta ";
            std::ostringstream ss;    ss << beta;   msg += ss.str();
            msg += "    ";
            output->toplog(msg);
        }
        
        
    }
    
    return beta;
    
}




// ---------------------------------------------------------------
// update the solution by executing the reactions involved in the fix
void Solution::update(int apos, double pV, int EVtype)
{
    
    if (msk->wplog){
        std::string msg;
        msg = "\tUPDATING SOLUTION at tempo ";
        std::ostringstream ss;
        ss << msk->tempo;   msg += ss.str()+"\n";    ss.str(""); ss.clear();
        msg += "Arguments are:  apos ="; ss << apos;   msg += ss.str();    ss.str(""); ss.clear();
        msg += "        pV =";             ss << pV;   msg += ss.str();    ss.str(""); ss.clear();
        output->toplog(msg);
    }
    
    // Compute number of chains, depending on mechanism for the fix
    
    int mid = fix->afKMCmid[apos];    // mechanism ID
    
    double nrv = 0;   // number of chains in particle.  NB: a double this fractions of molecules can be passed. This ensures correct volume changes associated to reaction and particle packing
    int chID = -1;  // chain ID: will be > -1 only if mechanism calls indeed a chain
    int rxid = -1;  // reaction ID: will be used now but later in loop too
    int nrxCH = -1; // number of reactions in series forming the chain
    
    if (strcmp((chem->mechstyle[mid]).c_str(),"allpar")==0 || strcmp((chem->mechstyle[mid]).c_str(),"allser")==0) {
        // compute number of chains depending on mechanism (allpar and allser should be the same, but shvab may be different due to radial growth idea
        if (chem->mechchain[mid]) {  //if reaction is a chain
            chID = chem->mechrcID[mid];
            nrv =  fabs(pV / chem->ch_dV_fgd[chID]);   //  made positive for either dissolution or nucleation events, i.e. irresepective of dV_fgd being negative or positive
            nrxCH = (chem->ch_rxID[chID]).size();
        }
        else{ //if instead it is a single reaction
            rxid = chem->mechrcID[mid];
            nrv = fabs(pV / chem -> rx_dV_fgd[rxid]);
            nrxCH = 1;
        }
    }
    else if (strcmp((chem->mechstyle[mid]).c_str(),"other_mechanisms_to_implement_here")==0) {
        
    }
    

    
    if (msk->wplog){
        std::string msg;
        std::ostringstream ss;
        msg = "Chain ID ="; ss << chID;   msg += ss.str()+"\t";    ss.str(""); ss.clear();
        msg += "nrv ="; ss << nrv;   msg += ss.str()+"\t";    ss.str(""); ss.clear();
        msg += "nrxCH =";             ss << nrxCH;   msg += ss.str()+"\n";    ss.str(""); ss.clear();
        output->toplog(msg);
    }
    
    
    if (EVtype==0) { //dissoltuion event
        SolidV -= pV;   // ATTENTION: ONLY TRUE FOR SPHERICAL PARTICLES!
    }
    else if (EVtype==1) { //nucleation event
        SolidV += pV;
    }
    PackF = SolidV/BoxV;
    
    
    if (strcmp(fix->afKMCsoutUL[apos].c_str(),"uniform")==0) {

        for (int i=0; i<nrxCH; i++) {
            
            int nrxpar = 1;     // number of reactions to be run in parallel per chain step (1 if not a chain)
            if (chem->mechchain[mid]) {
                rxid = chem->ch_rxID[chID][i];
                nrxpar = chem->ch_nrx[chID][i];
            }
            
            int nmols = chem->bkg_molID[rxid].size();  // number of molecule types involved in the background solution (to be created or removed depending on deletion or nucleation event type)
            
            double Vfrac = 1.;
            double dVfrac = 0.;  //to distribute molecules between box and dV. Default values assume afKMCsoutbox[apos] == box
            if (strcmp(fix->afKMCsoutbox[apos].c_str(),"box+dV")==0){
                Vfrac = SVol/(SVol+dVSVol);
                dVfrac = dVSVol/(SVol+dVSVol);
            }
            else if (strcmp(fix->afKMCsoutbox[apos].c_str(),"box")!=0){
                std::string msg = "ERROR: solution out should be \"box\" or \"box+dV\", instead "+fix->afKMCsoutbox[apos]+" was found for "+fix->afKMCname[apos];
                error->errsimple(msg);
            }
                
            // add/remove all types of bkg molecules altered by the reaction
            
            if (msk->wplog){
                std::string msg;
                msg = "Number of molecules in solution of:\n";
                output->toplog(msg);
            }
            
            for (int j=0; j<nmols; j++) {
                int molID = chem->bkg_molID[rxid][j];
                chem->mol_nins[molID] += (nrv * (double)nrxpar * chem->bkg_nmol[rxid][j] * Vfrac);
                chem->mol_nindV[molID] += (nrv * (double)nrxpar * chem->bkg_nmol[rxid][j] * dVfrac);
                
                if (msk->wplog){
                    std::string msg;
                    std::ostringstream ss;
                    msg = "- "; ss << molID;   msg += ss.str()+"\t";    ss.str(""); ss.clear();
                    msg += " --> in sol: "; ss << chem->mol_nins[molID] ;   msg += ss.str()+"\t";    ss.str(""); ss.clear();
                    msg += " --> in dV: "; ss << chem->mol_nindV[molID] ;   msg += ss.str()+"\n";    ss.str(""); ss.clear();
                    output->toplog(msg);
                }
            }
            // change of solution volume, concentrations, and volume of voids
            SVol += chem->rx_dV_bkg[rxid] * nrv * (double)nrxpar*Vfrac;
            dVSVol += chem->rx_dV_bkg[rxid] * nrv * (double)nrxpar*dVfrac;
            for (int j=0; j<nmols; j++) {
                int molID = chem->bkg_molID[rxid][j];
                chem->mol_cins[molID] = chem->mol_nins[molID]/SVol/nAvo/unitC;
                chem->mol_cindV[molID] = chem->mol_nindV[molID]/dVSVol/nAvo/unitC;
            }
        }
        
        // update voids in box and dV
        voidV = BoxV - SolidV - SVol;
        dVvoidV = dV - dVSVol;
    } else if (strcmp(fix->afKMCsoutUL[apos].c_str(),"local")==0) {
        // to be implemented: in this case only the subcommunnicator with managing the solution diffusion algorithm should do operations here. In this case, when compbeta is called by the fixes, it is again the subcomm manageing the solution that should be invoked. All to be implemented, so ERROR GIVEN FOR NOW
        std::string msg = "ERROR: local solution out for fix "+fix->afKMCname[apos]+" not implemented yet";
        error->errsimple(msg);
    } else {
        std::string msg = "ERROR: solution out should be \"uniform\", instead "+fix->afKMCsoutUL[apos]+" was found for "+fix->afKMCname[apos];
        error->errsimple(msg);
    }
}

// For now, only used by fix_nufeb.cpp
void Solution::updateconc(int apos, const std::vector<double>& dconc)
{
  double Vfrac = 1.;
  double dVfrac = 0.;  // to distribute molecules between box and dV. Default values assume afKMCsoutbox[apos] == box
  if (strcmp(fix->aCsoutbox[apos].c_str(),"box+dV")==0){
    Vfrac = SVol/(SVol+dVSVol);
    dVfrac = dVSVol/(SVol+dVSVol);
  }
  else if (strcmp(fix->aCsoutbox[apos].c_str(),"box")!=0){
    std::string msg = "ERROR: solution out should be \"box\" or \"box+dV\", instead "+fix->afKMCsoutbox[apos]+" was found for "+fix->afKMCname[apos];
    error->errsimple(msg);
  }
  for (int i=0; i<chem->Nmol; i++) {
    double dn = dconc[i] * SVol * nAvo * unitC; // jump in number of molecules based on jump in concentration
    double dnins = dn * Vfrac;
    double dnindV = dn * dVfrac;
    chem->mol_nins[i] += dnins;
    chem->mol_nindV[i] += dnindV;
    SVol += dnins * chem->mol_vapp[i];
    dVSVol += dnindV * chem->mol_vapp[i];
  }
  for (int i=0; i<chem->Nmol; i++) {
    chem->mol_cins[i] = chem->mol_nins[i]/SVol/nAvo/unitC;
    chem->mol_cindV[i] = chem->mol_nindV[i]/dVSVol/nAvo/unitC;
  }
  // update voids in box and dV
  voidV = BoxV - SolidV - SVol;
  dVvoidV = dV - dVSVol;
}

// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Solution::printall()
{
	fprintf(screen,"\n---------ALL ABOUT SOLUTION----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
