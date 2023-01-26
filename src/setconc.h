#ifndef SETCONC_H
#define SETCONC_H


#include "mpi.h"
#include "pointers.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>


#define MASTER 0
#define nAvo 6.022e23

namespace MASKE_NS {
    class Setconc : protected Pointers {
        
    public:
        
        Setconc(class MASKE *);
        ~Setconc();
        
        bool flag_setconc;   // true if setconc fix is active
        //std::vector<std::string> setnames;   //names of all user-defined setconc commands
        std::vector<std::string> molnames;   //names of molecule set by each setconc command
        std::vector<double> molconcs;       // concentrations of molecule to be set
        int vevery;       // integer with frequency of calls of commands
        bool ctr_flags;        // flags if counterions to be added in each command
        std::string ctr_mols;  // names of counterions
        std::string boxdV;  // string saying whether ions are to be fixed in the box only or also in dV
            
        std::vector<int> molID;   //IDs of molecule to be set
        int ctrID;     // IDs of counterion to be set
        
        void add_conc(std::string ,double);  // add set concentration
        void exec(void);   // sets the concentrations
        void printall();
        
    private:
        int me;
    };
}

#endif
