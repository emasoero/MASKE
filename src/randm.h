#ifndef RANDM_H
#define RANDM_H

#include "pointers.h"
#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include <unistd.h>   //just for the sleep() function
//#include <string>

#define MASTER 0

namespace MASKE_NS {
    class Randm : protected Pointers {
		
	public:
        
		Randm(class MASKE *);
		~Randm();
		
        int new_seme;   //random seed to be passed around among processor during seeding
        MPI_Status status;
        
        double pick01();
        void seedit(int);


        
        void printall();
        
    private:
        int me;
        int nprocs;     // number of processor in the universe (Comm_size)

	};
}

#endif
