#ifndef SPEC_H
#define SPEC_H

#include "pointers.h"
#include <IPhreeqc.hpp>
#include <string>
#include <vector>

#define MASTER 0

namespace MASKE_NS {
  class Spec : protected Pointers {
        
  public:
        
    Spec(class MASKE *);
    ~Spec();

    std::vector<std::string> specID; // IDs of Spec fixes in current subcomm
    std::vector<std::string> spec_type; // Type of speciation, currently only PHREEQC supported
    std::vector<int> spec_every; // frequency of each Spec in current subcomm
    std::vector<std::string> spec_database; // vector of databases
    std::vector<int> spec_solvents; // indices of solvent molecules
    std::vector<double> spec_masses; // molar masses in g/mol
    void add_spec(std::string,std::string,int,std::string,const std::vector<std::string>&,const std::vector<double>&);
    void dospec(int); // runs speciation
    void printall();

  private:
    IPhreeqc *IPhreeqc_ptr;
    int ipH, ipOH; // indices for H+ and OH- molecules
    std::map<std::string, double> residuals; // sum of speciation concentrations not considered in maske in mol 
  };
}

#endif
