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

namespace MASKE_NS {
    class Setconc : protected Pointers {
        
    public:
        
        Setconc(class MASKE *);
        ~Setconc();
        
        std::vector<std::string> setnames;   //names of all user-defined setconc commands
        std::vector<int> vevery;       // vector with frequency of calls of commands
        std::vector<std::string> molnames;   //names of molecule set by each setconc command
        std::vector<double> molconcs;       // concentrations of molecule to be set
        std::vector<bool> ctr_flags;        // flags if counterions to be added in each command
        std::vector<std::string> ctr_mols;  // names of counterions
        std::vector<std::string> vec_boxdV;  // vector saying whether ions are to be fixed in the box only or also in dV
            
        std::vector<int> molID;   //IDs of molecule to be set
        std::vector<int> ctrID;     // IDs of counterion to be set
        
        void add_conc(std::string , int, std::string ,double ,bool ,std::string ,std::string);  // add set concentration
        void exec(int);   // sets the concentrations
        void printall();
        
    private:
        int me;
    };
}

#endif
