#ifndef DT_NUCLEATE_H
#define DT_NUCLEATE_H

#include "pointers.h"
#include "mpi.h"
#include <string>
#include <vector>
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#define MASTER 0

namespace MASKE_NS {
	
	class DTnucleate : protected Pointers {
	public:
		//int narg;                    // # of command args
		//char **arg;                  // parsed args for command
		
		 DTnucleate(class MASKE *);
		~DTnucleate();
		
		bool isactive(int);
		void run(int);
        void printall();
        void addtrans(std::string,std::string,std::string,double,double,double);    // all processors will initially know all the nucleation events, that later will be allocated to subcommunicators via the allocSubcomm function
        void addtransSC(int);    // copies transition info into communicator-specific vectors
        
        int Ntrans;         // number of types of nucleation transitions defined by the user
        std::vector<std::string> Tnames;    // vector containing all the names of the transitions
        std::vector<std::string> scnames;   // vector containing all the names of the subcommunicators associated with each transition
        std::vector<double> freqeval;    // vector containing the frequenccy of transition evaluation, either in terms of time or steps
        std::vector<double> lasteval;    // vector containing the last step or time when the transition has been evaluated
        std::vector<double> starteval;    // vector containing the first step or time when the transition must be evaluated
        std::vector<double*> steptimeptr;   // vectors containing pointers to the current step or time of the system (to compare with the freuquency in the isactive() function)
        
        
        int Ntsc;         // number of nucleation transitions in current subcommunicator
       
        std::vector<std::string> TnamesSC;    // vector containing subcomm-specific transition names
        std::vector<double**> steptimeptrSC;   // vectors containing subcom-specific pointers to the global pointers to the to the current step or time of the system
        std::vector<double> lastevalSC;    // vector containing the subcomm-specific last step or time when the transition has been evaluated
        std::vector<double> startevalSC;    // vector containing the subcomm-specific first step or time when the transition must be evaluated
        std::vector<double> freqevalSC;    // vector containing the subcomm-specific frequenccy of transition evaluation, either in terms of time or steps
        std::vector<double> nextevalSC;    // vector containing the subcomm-specific step or time of next evaluation of transition, for each transition of this type

        
        
	private:
		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
	};
	
}

#endif
