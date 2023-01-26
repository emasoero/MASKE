#ifndef FIX_EMSTORELMP_H
#define FIX_EMSTORELMP_H

#include "pointers.h"
//#include <string>
#include <vector>
#include "mpi.h"
//#include <unistd.h>   //just for the sleep() function
//#include <string>
//#include <vector>
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#define MASTER 0
#define M_PI 3.14159265358979323846

namespace MASKE_NS {
	
	class Fix_EmstoreLMP : protected Pointers {
	public:
        Fix_EmstoreLMP(class MASKE *);
		~Fix_EmstoreLMP();
    
        
        // The quantities here are remembered and added when multiple fixes of same type are run on same subcomm
        //std::vector<int> IDsrt;         // vector storing the IDs of all delete events, EVEN IF PERTAINING TO DIFFERENT FIXES.   IDsrt is made of blocks, each containing sorted ID of particles invovled in a rej free KMC fix (any type) managed by the current subcomm.    Hence IDsrt will have size sum_j (Ng_j) where Ng_j is the number of particle involved in th j-th fix managed by the current subcomm.     IDs are sorted within each block (ie for each fix run by the subcomm), not overall.
        //std::vector<double> rates;         // vector storing the rate of all delete events, even if pertaining to different fixes

        // ----------------
        
        
        
        void execute(int);  // function executing fix of this type
        void printall();

        


        
        
	private:
        
        int me;     // id of the current processor (rank)
        
        /*MPI_Status status;

        // all the quantities below are "forgotten" after each fix-specific function is called.
        
        int natoms; // number of all atoms in lammps (certainly >= Ng)
        
        //----------------------------------------
        // quantities here are processor specific
        
        int key;    // rank of processor in its subcom. key = 0 will mean a submaster
        int mid;    // mechanism id to read mechanisms from fix.cpp
        double *aID;    // array with size natoms, with IDs from lammps but meaningful entries only for atoms in current processor
        double *aR, *aE;    // arrays with atom radius and energy, to be read from lammps
        std::vector<int> tID;   // proc-specific IDs
        std::vector<double> tR, tE;// same as aR and aE but with meaningful entries only for atoms in current processor
        int *aaID;      // array with size natoms, with IDs from lammps for ALL natoms
        int naP;      // number of atoms in current processor (naP <= Ng)
        int *tIDarr ;   // array of size naP with unsorted IDs of atoms in current processor only
        int *nID_each;     // array storing IDs in each processor: used by submaster for assembly and back
        double *rate_each;     // array storing rates in each processor

        
        //----------------------------------------
        // quantities here are assembled through all processors in subcomm, but are not kept going from one fix to another one of the same type on the same subcomm
        int Ng;     // number of atoms in group, viz particles involved and managed by the submaster for the current fix
        int nploc;  // number of processor in subcom
        int *IDuns;     // unsorted IDs of all Ng particles in fix. This is bloc-wise, with each block being the tIDarr of a processor
        double *unsRates;   // array in subcomm storing unsorted rates corresponding to IDuns above (size Ng)
        int *IDpos;  // array storing starting positions of local tID arrays in submaster's IDuns
        std::vector<int> UtoS;   //pointing each-proc unsorted IDuns to corresponding sorted IDsrt in global vector containing all rfreeKMC fixes in current subcomm.  Its size is Ng (number of particles in current fix).     Needed to send rates computed by each processor and assembled unsorted by submaster to global sorted-ID rate vector


        //----------------------------------------
        // internal functions follow
        void read_IDRE_lammps(int);    // read particle ID, Radius, and Energy from lammps, storing them an "aID", "aR", "aE" arrays with size natoms (all atoms in lammps) but meaningful entries only for particles in current processor
        void record_IDRE_proc(int);    // function to record particle ID, R, and E from aID.. to array tID.. arrays containing only the particles in the current processor
        void ids_to_submaster(int);     // assembles all IDs from each processor into block-wise array in submaster, where each block is the tIDarr in a processor
        void submaster_sort_IDs(int);   // submaster sorts IDuns and places them in the block corresponding to the current fix in the global block-wise-sorted vector IDsrt (see in public part of defs)
        void submaster_map_ID(int);    // the submaster maps unsorted IDuns to corresponding IDsrt entry in global blockwise-sorted vector containing all fixes run by current subcomm
        void comp_rates_allpar(int);    // each proc computes deletion rates
        void rates_to_submaster(int);   // submaster assembles rates in vector corresponding to unsorted IDs
        */
        
        
	};
	
}

#endif
