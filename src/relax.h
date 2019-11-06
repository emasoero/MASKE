#ifndef RELAX_H
#define RELAX_H


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
    class Relax : protected Pointers {
        
    public:
        
        Relax(class MASKE *);
        ~Relax();

        std::vector<std::string> rlxID; // IDs of Relax fixes in current subcomm
        std::vector<int> rlx_every; //frequency of each Relax in current subcomm
        std::vector<std::string> rlx_relaxers; // vector of relaxer modes (minimize, nvt, ..)
        std::vector<std::string> rlx_string; // string to run in minimize or nvt
        std::vector<std::string> rlx_style; // stlye of minimizers
        std::vector<std::string> rlx_modify; // string to run in min_modify
        
        void add_rlx(std::string,int,std::string,std::string,std::string,std::string);  // add relaxer
        void dorelax(int);   // runs the relaxer
        void printall();
        
    private:
        int me;
    };
}

#endif
