#ifndef MASKE_H
#define MASKE_H

#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "mpi.h"
#include <unistd.h>   //just for the sleep() function
#include <vector>



/*#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
*/

#define MASTER 0


namespace MASKE_NS {
    
    //using namespace LAMMPS_NS;
    
    
	class MASKE {
		
        
        
	public:
		
		class Memory *memory;          // memory allocation functions
        class Error *error;
        class Inputmsk *inputmsk;
        class Universe *universe;
        class DTnucleate *nucleate;
        class LammpsIO *lammpsIO;
        class Simbox *simbox;
        class Particles *particles;
        class Chemistry *chem;
        class Interactions *interact;
        class Solution *solution;
        class Fix *fix;
        class Krun *krun;
        class Randm *randm;
        class Fix_delete *fix_del;
        class Output *output;
        class Fix_Cfoo *fix_cfoo;
        class Relax *relax;
        class Fix_nucleate *fix_nucl;
        class Block *block;
        class Store *store;
		
        
        int Nsteps;     // the total number of fundametnal simulation steps to run (each transition is called every multiple of this step)
        int step;       // the current simulation step
        double doublestep;  // same as step, but in double version, as needed to calculate the frequency of invoking transition (see isasctive functions in transition files)
        double tempo;   // current simulation time
        double kB;      // Boltzmann constant (user-defined in input, otherwise initialised to 1)
        double hpl;      // Planck constant (user-defined in input, otherwise initialised to 1)
        
		FILE *screen;                  // screen output
        bool wplog;      // if true, each processor writes a processor specific log for debug
        FILE *plog;                  // processor-specific log file output
        std::string plogfname;
        bool wthermo;      // if true, the MASTER writes the thermo file
        FILE *thermo;                  // thermo outup file: written by master only
        std::string th_fname;

        std::vector<int> Rtypes;    // vector of types associated to real particles
        std::vector<int> Ttypes;    // vector of types associated to trial particles
        
		MASKE(int,char**);  // constructor
		~MASKE();           // destructor
		
		void printall();
        void mainloop();   //OBSOLETE
        
    private:
        int me;
	};
}

#endif
