#ifndef SOLUTION_H
#define SOLUTION_H

#include "pointers.h"
#include <string>
#include <stdlib.h>
#include <vector>
#include <math.h>       /* for the pow and sqrt and fabs functions used here
#include <unistd.h>   //just for the sleep() function

/*#include "mpi.h"

#include <stdio.h>
#include <iostream>
#include <fstream>


#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
*/


#define MASTER 0
#define M_PI 3.14159265358979323846
#define nAvo 6.022e23

namespace MASKE_NS {
    
	class Solution : protected Pointers {
	public:
		
		 Solution(class MASKE *);
		~Solution();
		
        std::vector<std::string> molins;    // vector with types of molecules in solution
        std::vector<double> conc, concdV;           // vector with concentration of molecules in solution in box and in dV (when aniso dV implemented, a 2D vector should be used instead of concdV)
        std::vector<double> nmol, nmoldV;           // vector with number of molecules in solution in box and in dV (when aniso dV implemented, a 2D vector should be used instead of nmoldV)
        std::vector<double> molID;         // vector with ids of molecules in solution (to point to molecule in chemistry.cpp vectors)
        std::string soltype;                // type of solution: uniform vs lattice etc..
        double Temp;                        // Temperature: for now can be only uniform
        std::string dVtype;     // iso or aniso
        double voidV, dVvoidV ;         // volume of voids in box and in dV (when aniso dV implemented, a vector should be used instead of dVvoidV)

        double DH_A, DH_B;                  // Debye-Huckel constants A and B, depending on the solvent. Must be user-defined in input.dat
        
        double BoxV, dV, PackF, SVol, dVSVol,SolidV;  // box volume,  additional reservoir volume of solution (when aniso dV implemented, a vector should be used instead), solid fraction, and volume of solution in box and in dV (when aniso dV implemented, a vector should be used instead of dVvoidV)
        double unitC;   // conversion factor to go from user's unit volume (e.g. nm3) to litre
        
        MPI_Status status;
        
        void addmol(std::string  , double );
        void printall();
        void computeNmol();     // function to compute the number of molecules in solution in box and in dV
        double compbeta(int,bool,std::string,std::string);
        void update(int,double,int);
	void updateconc(int,const std::vector<double>&);
        
        double compQ(int,std::string,std::string,std::string);   // computes activity products of various things depending on second input argument. (1) If it is "reac", these are background reactants in a precipitation or dissolution reaction, i.e. only  negative bkg terms in the reaction; (2) if "prod", products in a precipitation or dissolution reaction, i.e. only positive bkg terms in the reaction.    NOTICE that fgd terms are neglected, which is problematic if some charged molecules in solution are represented as fgd particles, as they would not be accounted for in the ionic strength.
        
        int me;     // id of the current processor (rank)

        
	private:

        
	};
	
}

#endif
