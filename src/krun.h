#ifndef KRUN_H
#define KRUN_H

#include "pointers.h"
//#include <string>
#include <vector>
#include "mpi.h"
#include "stdlib.h"
#include <unistd.h>   //just for the sleep() function
#include <math.h>       /* for the log (natural log) function */
//#include <string>
//#include <vector>
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#define MASTER 0

namespace MASKE_NS {
	
	class Krun : protected Pointers {
	public:
        Krun(class MASKE *);
		~Krun();
        
        bool reset_event_vec;
        std::vector<double> RFKMCrates;  // current rates of all RF-KMC event in current subproc
        std::vector<double> RFKMCratesIN;  // rates of all RF-KMC event in current subproc integrated over time until now
        //std::vector<double> CumRFKMCratesIN;  // cumulative of RFKMCratesIN, needed only when choosing the RF-KMC event to perform when the time of occurrence is reached or passed
        bool fMC_justreset; //records when the vectors of rej-free events have been cleared. Used by process-specific fix_type.cpp to choose when to assemble new ID, energy, rate, etc vectors
        
        double QkTint;   // integral of total cum rate of all RF-KMc events from time of last event (t0) until current time (tempo)
        bool resetQkTint;   // if true, QkTint is reset to zero before executing krun. This may be useful when the user introduces new KMC fixes from one krun to the next
        
        void proceed(double);
        void printall();
        
        
	private:
		//std::string fname;              // open inputcprs file
        int me;     // id of the current processor (rank)
        MPI_Status status;
        MPI_Request request;

	std::vector<double> nufeb_buf; // buffer for exchange of atoms when running nufeb
	};
	
}

#endif
