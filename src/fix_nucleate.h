#ifndef FIX_NUCLEATE_H
#define FIX_NUCLEATE_H

#include "pointers.h"
//#include <string>
#include <vector>
#include "mpi.h"
#include <unistd.h>   //just for the sleep() function
#include <lmptype.h>

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
	
	class Fix_nucleate : protected Pointers {
	public:
        Fix_nucleate(class MASKE *);
		~Fix_nucleate();
    
        
        // The quantities here are remembered and added when multiple fixes of same type are run on same subcomm
       
         std::vector<int> IDsrt;         // vector storing the IDs of all delete events, EVEN IF PERTAINING TO DIFFERENT FIXES.   IDsrt is made of blocks, each containing sorted ID of particles invovled in a rej free KMC fix (any type) managed by the current subcomm.    Hence IDsrt will have size sum_j (Ng_j) where Ng_j is the number of particle involved in th j-th fix managed by the current subcomm.     IDs are sorted within each block (ie for each fix run by the subcomm), not overall.
        std::vector<double> rates;         // vector storing the rate of all nucleation events, even if pertaining to different fixes
        std::vector<double> ratesT;         // vector storing the time-integrated  rates of all nucleate events, even if pertaining to different fixes
        std::vector<double> CratesT;         // vector cumulating ratesT
        std::vector<int>cumNPfix,fixaID;    //vectors with cumulated number of particles per nucleation-type fix, and absolute ID of that fix
        
        int EVpID;  //  particle ID for the chosen nucleate event
        int EVpTYPE;   //  particle type for the chosen nucleate event
        int EVafixID; // absolute fix ID of chosen particle for nucleate event
        double EVpDIAM, EVpXpos, EVpYpos, EVpZpos;  // diameter and coordinates of chosen particle for nucleate event
         double partVol;     // particle volume to be communicated around when an event is chosen and executed. Needed then to update the solution after event execution

        // ----------------
        
        
        
        void init(int);       // initialise the fix: here creates trial particles on lattice in region and counts them
        void delete_trial(int);       // delete the trial particles in a nucleate fix
	void reset_trial_pos(int); // reset position of trial particles
        void sample(int);       // function computing the positions, IDs, and rates of all events and adding them to the vectors above
        double CrateAll();      // function computing cumulated unsorted rate vectors and returning the cumulated rate value
        double CratesTupdate(double);   // //updates the vector of time integrals of  cumulated rates, CratesT
        void clearvecs();
        
        void findEV(double);    // find particle to be nucleated based on cumulative rate extracted randomly. Also extracts particle position and shape parameters to be then accessed from Krun.cpp
        void extractPROPS(int,int);    // extract properties attached to particle to be nucleated: position and shape. Arguemnts are particle id and absolute fix id
        void execute(int, int , int, double ,double , double , double);    // execute nucleation event of particle with given ID, all-fix ID, type, diameter and X,Y,Z coordinates

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
        
        
         
        double *aID;    // array with size natoms, with IDs from lammps but meaningful entries only for atoms in current processor
        //double *aR,
        double *atype;
        double *aE;    // arrays with atom radius and energy, to be read from lammps
        double *aR;    // arrays with atom radius and energy, to be read from lammps
        std::vector<int> tID;   // proc-specific IDs
        std::vector<double> tR; // same as aR and aE but with meaningful entries only for atoms in current processor
        std::vector<double> tE;// same as aR and aE but with meaningful entries only for atoms in current processor
        double *tCF; // coverage fractions of atoms in current processor
        double *tGM; // maximum interfacial energy change of atoms in current processor
        //int *aaID;      // array with size natoms, with IDs from lammps for ALL natoms
        int naP;      // number of atoms in current processor (naP <= Ng)
        int *tIDarr ;   // array of size naP with unsorted IDs of atoms in current processor only
        int *nID_each;     // array storing IDs in each processor: used by submaster for assembly and back
        double *rate_each;     // array storing rates in each processor
        double **locLMP;   // array to extract content of compute property/local from LAMMPS, used in "micro" style to obtain IDs and types of interacting pairs of atoms from neighbour list
        double *aDIST;   // vector to extract content of compute pair/local from LAMMPS, used in "micro" style to obtain distances between interacting pairs of atoms from neighbour list
        int nlocR; // number of rows from LAMMPS property/local (neighbour list)
        int *nlocR_each;     // array storing numcber of local array rows in each processor: used by submaster for assembly and back (only used for "micro" style)
        
        // Block of variables below ar for coveraage fraction calculations in micro mechanism
        //********
        int *IDar;   // IDs of all real particles
        double *Rar;   // Radii of all real particles
        std::vector<int> tIDar;   // temporary vector with ID of all real particles in current proc
        std::vector<double> tRar;// temporary vector with Radii of all real particles in current proc
        int nR;      // number of real atoms in current fix in current processor (for coveraage fraction calculations in micro mechanism)
        int *nIDar_each;     // array storing number of real particle IDs in each processor
        //********
        
        
        //----------------------------------------
        // quantities here are assembled through all processors in subcomm, but are not kept going from one fix to another one of the same type on the same subcomm
        int Ng;     // number of atoms in group, viz particles involved and managed by the submaster for the current fix
        int nploc;  // number of processor in subcom
        int *IDuns;     // unsorted IDs of all Ng particles in fix. This is bloc-wise, with each block being the tIDarr of a processor
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

        // Block of variables below ar for coveraage fraction calculations in micro mechanism
        //********
        int *IDarpos;  // array storing starting positions of local IDar arrays in submaster's IDaruns
        int *IDaruns;     // unsorted IDs of all real particles. This is bloc-wise, with each block being the IDar of a processor
        double *Raruns;     // unsorted Radii of all real particles. This is bloc-wise, with each block being the IDar of a processor
        int totRp;  // total number of real particles in all processors in LAMMPS
        std::vector<int> IDarsrt;     // sorted IDs of all real particles, only in submaster
        std::vector<double> Rarsrt;     // sorted Radii of all real particles, only in submaster
        //********

        int natoms0; // number of all atoms in lammps when initial position of trial particles is recorded
        int *id0;    // id of all LAMMPS atoms (real and trial) in this subcomm, when init() is called (used to later reset the initial positions of trial particles)
        int *typ0; // initial types
        double *xyz0; // initial positions
        
	/*
     int nmax;
	int nlocal0;
	LAMMPS_NS::tagint *tag0; // initial tags
	int *type0; // initial types
	double **x0; // initial positions
     */
	 
        //----------------------------------------
        // internal functions follow
        /*
        void read_IDRE_lammps(int);    // read particle ID, Radius, and Energy from lammps, storing them an "aID", "aR", "aE" arrays with size natoms (all atoms in lammps) but meaningful entries only for particles in current processor
        void record_IDRE_proc(int);    // function to record particle ID, R, and E from aID.. to array tID.. arrays containing only the particles in the current processor
        */
        void ids_to_submaster(int);     // assembles all IDs from each processor into block-wise array in submaster, where each block is the tIDarr in a processor
        void pair_arr_to_submaster(int);     // each proc communicates its portion of pair array (from neighbour list) to the submaster - used by style "micro" only, to compute coverage areas

        void submaster_comp_cover(int);     // The submaster computes coverage areas of particles in the current fix and communicates them back to the other processors in its subcomm
        void cover_from_submaster(int);     // each processor in the subcommunicator receives its chunk of unsorted coverage area fractions from the submaster
        
        void submaster_sort_IDs(int);   // submaster sorts IDuns and places them in the block corresponding to the current fix in the global block-wise-sorted vector IDsrt (see in public part of defs)
        
        void submaster_map_ID(int);    // the submaster maps unsorted IDuns to corresponding IDsrt entry in global blockwise-sorted vector containing all fixes run by current subcomm
        
        void comp_rates_allpar(int);    // each proc computes deletion rates
        void comp_rates_micro(int);    // each proc computes deletion rates
        
         void rates_to_submaster(int);   // submaster assembles rates in vector corresponding to unsorted IDs
        /*
         void comp_rates_allser(int);    // each proc computes deletion rates for an all-series mechanism

        */
        
        
	};
	
}

#endif
