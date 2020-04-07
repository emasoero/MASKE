#ifndef FIX_NUFEB_H
#define FIX_NUFEB_H

#include "pointers.h"
#include "mpi.h"

#include <fix_bio_kinetics.h>

#include <vector>

#define MASTER 0
#define M_PI 3.14159265358979323846

namespace MASKE_NS {
	
  class Fix_nufeb : protected Pointers {
  public:
    Fix_nufeb(class MASKE *);
    ~Fix_nufeb();
    
    void init(int pos);
    void setup(int);
    void setup_exchange(int);
    double getDT(int);       // function computing time increment for the continuum process
    void printall();
    void execute(int,int);
    void exchange(int,int);
    
    int group; // bacteria group

    int ncells;
    int first_flag;
    int setup_flag;
    int setup_exchange_flag;
    std::vector<double> sublo, subhi;
    std::vector<bool> intersect;
    std::vector<int> procs;
    std::vector<int> nsend, nrecv;
    std::vector<double> buf;

    LAMMPS_NS::FixKinetics *kinetics;
  };
}

#endif
