#ifndef FIX_DELETE_H
#define FIX_DELETE_H

#include "pointers.h"
//#include <string>
#include <vector>
#include "mpi.h"
#include <unistd.h>   //just for the sleep() function

//#include <string>
//#include <vector>
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#define MASTER 0
#define M_PI 3.14159265358979323846
#define nAvo 6.022e23

namespace MASKE_NS {
	
	class Fix_delete : protected Pointers {
	public:
        Fix_delete(class MASKE *);
		~Fix_delete();
    
        
        // The quantities here are remembered and added when multiple fixes of same type are run on same subcomm
        std::vector<int> IDsrt;         // vector storing the IDs of all delete events, EVEN IF PERTAINING TO DIFFERENT FIXES.   IDsrt is made of blocks, each containing sorted ID of particles invovled in a rej free KMC fix (any type) managed by the current subcomm.    Hence IDsrt will have size sum_j (Ng_j) where Ng_j is the number of particle involved in th j-th fix managed by the current subcomm.     IDs are sorted within each block (ie for each fix run by the subcomm), not overall.
        std::vector<double> rates;         // vector storing the rate of all delete events, even if pertaining to different fixes
        std::vector<double> ratesT;         // vector storing the time-integrated  rates of all delete events, even if pertaining to different fixes
        std::vector<double> CratesT;         // vector cumulating ratesT
        std::vector<int>cumNPfix,fixaID;    //vectors with cumulated number of particles per delete-type fix, and absolute ID of that fix

        int EVpID;  //  particle ID for the chosen delete event
        int EVafixID; // absolute fix ID of chosen particle for delete event
        double partVol;     // particle volume to be communicated around when an event is chosen and executed. Needed then to update the solution after event execution

        // ----------------
        
        
        void init(int);       // function initialising the fix: here counts the number of particles to be deleted by current fix
        void sample(int);       // function computing the IDs and rates of all events and adding them to the vectors above
        double CrateAll();      // function computing cumulated unsorted rate vectors and returning the cumulated rate value
        double CratesTupdate(double);   // //updates the vector of time integrals of  cumulated rates, CratesT
        void clearvecs();
        void findEV(double);    // find particle to be deleted based on cumulative rate extracted randomly
        void execute(int);       // function deleting a particle in lammps with ID passed as argument
        void add_Ng_apos(int,int);  // function adding elements to cumulative vector of number of particles in fix, and associated absolute ID of fix
        
        
        void printall();

    


	private:
        MPI_Status status;

        // all the quantities below are "forgotten" after each fix-specific function is called.
        
        int natoms; // number of all atoms in lammps (certainly >= Ng)
        int nlocal; // number of atoms in current processor
        
        //----------------------------------------
        // quantities here are processor specific
        int me;     // id of the current processor (rank)
        int key;    // rank of processor in its subcom. key = 0 will mean a submaster
        int mid;    // mechanism id to read mechanisms from fix.cpp
        double *aR, *aE;    // arrays with atom radius, energy, and type to be read from lammps
        double *atype;
        double *aID;    // array with size natoms, with IDs from lammps but meaningful entries only for atoms in current processor
        std::vector<int> tID;   // proc-specific IDs
        std::vector<double> tR, tE;// same as aR and aE but with meaningful entries only for atoms in current processor
        double *tCF; // coverage fractions of atoms in current processor
        double *tGM; // maximum interfacial energy change of atoms in current processor
        int naP;      // number of atoms in current fix in current processor
        int *tIDarr ;   // array of size naP with unsorted IDs of atoms in current processor only
        int *nID_each;     // array storing number of IDs in each processor: used by submaster for assembly and back
        double *rate_each;     // array storing rates in each processor
        double **locLMP;   // array to extract content of compute property/local from LAMMPS, used in "micro" style to obtain IDs and types of interacting pairs of atoms from neighbour list
        double *aDIST;   // vector to extract content of compute pair/local from LAMMPS, used in "micro" style to obtain distances between interacting pairs of atoms from neighbour list
        int nlocR; // number of rows from LAMMPS property/local (neighbour list)
        int *nlocR_each;     // array storing numcber of local array rows in each processor: used by submaster for assembly and back (only used for "micro" style)

        
        //----------------------------------------
        // quantities here are assembled through all processors in subcomm, but are not kept going from one fix to another one of the same type on the same subcomm
        int Ng;     // number of atoms in group, viz particles involved and managed by the submaster for the current fix
        int nploc;  // number of processor in subcom
        int *IDuns;     // unsorted IDs of all Ng particles in fix. This is bloc-wise, with each block being the tIDarr of a processor
        double *Runs;    // unsorted radii of all Ng particles in fix
        double *CFuns;    // coverage fraction area of Ng particles in the fix
        double *GMuns;    // maximum interfacial enrgy change of Ng particles in the fix
        bool *fGMuns;    // is this the first interfacial enrgy change computed for each of the Ng particles in the fix?
        double *unsRates;   // array in subcomm storing unsorted rates corresponding to IDuns above (size Ng)
        int *IDpos;  // array storing starting positions of local tID arrays in submaster's IDuns
        std::vector<int> UtoS;   //pointing each-proc unsorted IDuns to corresponding sorted IDsrt in global vector containing all rfreeKMC fixes in current subcomm.  Its size is Ng (number of particles in current fix).     Needed to send rates computed by each processor and assembled unsorted by submaster to global sorted-ID rate vector
        std::vector<int> StoU;   //pointing sorted IDs in current fix_delete (pos) to position in unsorted list assembled from processors on submaster. NOTE: the IDsrt vector contains all IDs for all fix_deletes on this subcomm. Instead , the StoU map only refers to IDs in current fix delete
        int *SARpos;  // "Submaster ARray position" vector, storing starting positions of local arrays (from neighbour list) in submaster's ARuns and Duns
        double **SAR;  //Submaster ARray, bloc-wise with each block being the locLMP of a processor
        double *Dsub;    // vector on the submaster gathering all interaction distances from lammps pair/local
    



        //----------------------------------------
        // internal functions follow
        void ids_to_submaster(int);     // assembles all IDs from each processor into block-wise array in submaster, where each block is the tIDarr in a processor
        void submaster_sort_IDs(int);   // submaster sorts IDuns and places them in the block corresponding to the current fix in the global block-wise-sorted vector IDsrt (see in public part of defs)
        void submaster_map_ID(int);    // the submaster maps unsorted IDuns to corresponding IDsrt entry in global blockwise-sorted vector containing all fixes run by current subcomm
        void comp_rates_allpar(int);    // each proc computes deletion rates
        void rates_to_submaster(int);   // submaster assembles rates in vector corresponding to unsorted IDs
        void comp_rates_allser(int);    // each proc computes deletion rates for an all-series mechanism
        void comp_rates_micro(int);     // each proc computes deletion rates for a "micro" type of mechanism
        void pair_arr_to_submaster(int);     // each proc communicates its portion of pair array (from neighbour list) to the submaster - used by style "micro" only, to compute coverage areas
        void submaster_comp_cover(int);     // The submaster computes coverage areas of particles in the current fix and communicates them back to the other processors in its subcomm
        void cover_from_submaster(int);     // each processor in the subcommunicator receives its chunk of unsorted coverage area fractions from the submaster
        
        
	};
	
}

#endif
