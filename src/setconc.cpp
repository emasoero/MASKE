#include "setconc.h"
#include "chemistry.h"
//#include "universe.h"
//#include "lammpsIO.h"
#include "error.h"
#include "solution.h"

#include <string.h>
#include "memory.h"

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
void Setconc::add_conc(std::string mname,double mconc)
{   
    
    molnames.push_back(mname);
    molconcs.push_back(mconc);
    // vevery.push_back(every);
    // ctr_flags.push_back(flag_ctr);
    // ctr_mols.push_back(cmol);
    // vec_boxdV.push_back(boxdV);
    
    // find ID corresponding to names of molecule to be set
    molID.push_back(-1);
    // ctrID.push_back(-1);
    for (int i=0; i<chem->molnames.size(); i++){
        if( strcmp(mname.c_str(),chem->molnames[i].c_str())==0) molID[molID.size()-1] = i;
        // if( setconc->ctr_flags && strcmp(setconc->ctr_mols,chem->molnames[i].c_str())==0) setconc->ctrID = i;
    }
    if (molID.back()==-1){
        std::string msg = "ERROR: molecule "+mname+" requested in setconc command not found\n";
        error->errsimple(msg);
    }
    //if (setconc->ctr_flags && ctrID.back()==-1){
    //    std::string msg = "ERROR: counterion "+setconc->ctr_mols+" requested in setconc command not found\n";
    //    error->errsimple(msg);
    //}
    
    //  // check that defined ion and counterion have opposite charges (or that molecule has zero charge)
    // if (setconc->ctr_flags){
    //    error if molecule and counterion charge have same sign
    //    if (chem->mol_z[molID.back()] * chem->mol_z[setconc->ctrID] > 0.){
    //        std::string msg = "ERROR: in setconc, counterion "+setconc->ctr_mols+" has same charge as molecule "+mname+" \n";
    //        error->errsimple(msg);
    //    }
    
    //  // error if molecule is charged and counterion is not
    //    if (chem->mol_z[molID.back()] != 0  && chem->mol_z[setconc->ctrID]==0.){
    //        std::string msg = "ERROR: in setconc, counterion "+setconc->ctr_mols+" is neutral, but molecule "+mname+" is charged\n";
    //        error->errsimple(msg);
    //    }
    //}
}


// ---------------------------------------------------------------
// set concentration of given molecule and possibly counterion
 void Setconc::exec(void)
{   
    // fprintf(screen,"\n DEBUG: Initial SolVol = %e \n",solution->SVol);
    // fprintf(screen,"\n\n DEBUG setconc1: Proc %d, concentration of nutrient (%d) mapped to molecule (%d) is %e \n",me, chem->mol_nufeb[3], 3, chem->mol_cins[3]);
    // sleep(1);
    double temp =0;
        for(int i=0;i<chem->Nmol;i++){
        
        
            temp+=chem->mol_nins[i]*chem->mol_vapp[i];
        }
    // fprintf(screen,"\n DEBUG: Check SolVol = %e \n",temp);


    int up1;
    up1=molnames.size()+1;
//    double** A;
//    memory->create(A,up1,up1);

    double** A = new double*[up1];
    for(int i = 0; i < up1; i++){
        A[i] = new double[up1];
    }

    double* b;
    b = new double[up1];

    for (int i=0; i<up1-1; i++){
        for (int j=0; j<up1-1; j++){
            if (i==j) A[i][j]=1-molconcs[i]*chem->mol_vapp[molID[j]]* nAvo * solution->unitC;
            else A[i][j]=-molconcs[i]*chem->mol_vapp[molID[j]]* nAvo * solution->unitC;
        }
        A[i][up1-1]=-molconcs[i]*chem->mol_vapp[ctrID]* nAvo * solution->unitC;
    }


    for(int j=0;j<up1-1;j++) A[up1-1][j]=chem->mol_z[molID[j]];
    A[up1-1][up1-1]=chem->mol_z[ctrID];

    double sum_volknown=0;
    double sum_charge=0;


    for(int i=0;i<chem->Nmol;i++){
        
        bool skip=false;
        
        for(int j=0;j<molID.size();j++){
            if (i==molID[j]) skip=true;
        }
        if (i==ctrID) skip=true;
            
        if(skip==false) {
            sum_volknown+=chem->mol_nins[i]*chem->mol_vapp[i];
            sum_charge-=chem->mol_nins[i]*chem->mol_z[i];
        }

    }


    // fprintf(screen,"\n DEBUG: KnownVol pre-gauss = %e \n",sum_volknown);

    double sum_volunknown=0;
    for(int i=0; i<up1-1;i++){
        sum_volunknown+=chem->mol_nins[molID[i]]*chem->mol_vapp[molID[i]];
    }
    // fprintf(screen,"\n DEBUG: UnknownVol pre-gauss = %e \n",sum_volunknown);

    for(int i=0; i<up1-1;i++){

        b[i]=molconcs[i]*sum_volknown * nAvo * solution->unitC;
    }
    b[up1-1]=sum_charge;
  
    for (int j = 0; j<up1-1; j++){
        for (int i = j+1; i<up1; i++){
                double m = A[i][j]/A[j][j];
            for (int k=0;k<up1;k++) A[i][k] = A[i][k] - m*A[j][k];
            b[i] = b[i] - m*b[j];
        }
        
    }


    double* x; //array of unknown values
    x = new double[up1];
    for (int k=0;k<up1;k++) x[k] = 0; //initialise to zero

    x[up1-1] = b[up1-1]/A[up1-1][up1-1];               
    for (int i = up1-2; i>-1; i--){                   
        double  sum = 0;
        for (int j = up1-1; j>i; j--){                
            sum += A[i][j]*x[j];    
        }
        x[i] = (b[i]- sum)/A[i][i];
    }
    

    //Assign the unknown values to the corresponding number of ion values and concentration values in chemistry 
    for(int i=0; i<up1-1;i++){
        // chem->mol_cins[molID[i]] = molconcs[i];
        chem->mol_nins[molID[i]] = x[i];
    }
    chem->mol_nins[ctrID]=x[up1-1];
    
    //Calculate volume of solution


    sum_volunknown=0;
    for(int i=0; i<up1-1;i++){
        sum_volunknown+=chem->mol_nins[molID[i]]*chem->mol_vapp[molID[i]];
    }


    // fprintf(screen,"\n DEBUG: UnknownVol post-gauss = %e \n",sum_volunknown);
    // fprintf(screen,"\n DEBUG: KnownVol post-gauss = %e \n",sum_volknown);
    double sum_volctr;
    sum_volctr=chem->mol_nins[ctrID]*chem->mol_vapp[ctrID];
    // fprintf(screen,"\n DEBUG: Ctrvol post-gauss =  %e \n",sum_volctr);

    solution->SVol=sum_volknown+sum_volunknown+sum_volctr;
    // fprintf(screen,"\n DEBUG: SolutionVol post-gauss =  %e \n",solution->SVol);


    //Recalculate the concentrations of all species
    for(int i=0;i<chem->Nmol;i++){
        chem->mol_cins[i]=chem->mol_nins[i]/(solution->SVol * nAvo * solution->unitC);
    }

    //starting boxdV

    if (strcmp(boxdV.c_str(),"box+dV")==0) {
        
        sum_volknown=0;
        sum_charge=0;

        for(int i=0;i<chem->Nmol;i++){
            
            bool skip=false;
            
            for(int j=0;j<molID.size();j++){
                if (i==molID[j]) skip=true;
            }
            if (i==ctrID) skip=true;
                
            if(skip==false) {
                sum_volknown+=chem->mol_nindV[i]*chem->mol_vapp[i];
                sum_charge-=chem->mol_nindV[i]*chem->mol_z[i];
            }

        }
        for(int i=0; i<up1-1;i++){

            b[i]=molconcs[i]*sum_volknown * nAvo * solution->unitC;
        }
        b[up1-1]=sum_charge;

        
        for (int j = 0; j<up1-1; j++){
            for (int i = j+1; i<up1; i++){
                    double m = A[i][j]/A[j][j];
                for (int k=0;k<up1;k++) A[i][k] = A[i][k] - m*A[j][k];
                b[i] = b[i] - m*b[j];
            }
            
        }


        for (int k=0;k<up1;k++) x[k] = 0; //initialise to zero

        x[up1-1] = b[up1-1]/A[up1-1][up1-1];               
        for (int i = up1-2; i>-1; i--){                   
            double  sum = 0;
            for (int j = up1-1; j>i; j--){                
                sum += A[i][j]*x[j];    
            }
            x[i] = (b[i]- sum)/A[i][i];
        }

        //Assign the unknown values to the corresponding number of ion values and concentration values in chemistry 
        for(int i=0; i<up1-1;i++){
            // chem->mol_cins[molID[i]] = molconcs[i];
            chem->mol_nindV[molID[i]] = x[i];
        }
        chem->mol_nindV[ctrID]=x[up1-1];
        
        //Calculate volume of solution

        sum_volunknown=0;
        for(int i=0; i<up1-1;i++){
            sum_volunknown+=chem->mol_nindV[molID[i]]*chem->mol_vapp[molID[i]];
        }
        
        sum_volctr=chem->mol_nindV[ctrID]*chem->mol_vapp[ctrID];

        solution->SVol=sum_volknown+sum_volunknown+sum_volctr;

        //Recalculate the concentrations of all species
        for(int i=0;i<chem->Nmol;i++){
            chem->mol_cindV[i]=chem->mol_nindV[i]/(solution->SVol * nAvo * solution->unitC);
        }
    }  

    for (int i=0; i<up1; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;



    

/*   double c_old = chem->mol_cins[molID[i]];
     double cdV_old = chem->mol_cindV[molID[i]];
     double n_old = chem->mol_nins[molID[i]];
     double ndV_old = chem->mol_nindV[molID[i]];
     
     chem->mol_cins[molID[i]] = molconcs[i];
     chem->mol_nins[molID[i]] = chem->mol_cins[molID[i]] * solution->SVol * nAvo * solution->unitC;
     
     if (strcmp(vec_boxdV[i].c_str(),"box+dV")==0) {
         chem->mol_cindV[molID[i]] = molconcs[i];
         chem->mol_nindV[molID[i]] = chem->mol_cindV[molID[i]] * solution->dVSVol * nAvo * solution->unitC;
     }
     
     // same for counterion
     if (ctr_flags[i]){
         double dn;
         dn =  (chem->mol_nins[molID[i]] - n_old) * (-chem->mol_z[molID[i]]/chem->mol_z[ctrID[i]]);
         chem->mol_cins[ctrID[i]] += dn/solution->SVol/nAvo/solution->unitC;
         chem->mol_nins[ctrID[i]] += dn;
         
         if (strcmp(vec_boxdV[i].c_str(),"box+dV")==0) {
             dn =  (chem->mol_nindV[molID[i]] - ndV_old) * (-chem->mol_z[molID[i]]/chem->mol_z[ctrID[i]]);
             chem->mol_cindV[ctrID[i]] += dn/solution->dVSVol/nAvo/solution->unitC;
             chem->mol_nindV[ctrID[i]] += dn;
         }
     }*/
    
}


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Setconc::printall()
{   
    fprintf(screen,"\n---------ALL ABOUT SETCONC----------\n");
    //fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
    fprintf(screen,"---------------------------------------\n\n");
    
}
