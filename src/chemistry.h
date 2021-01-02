#ifndef chemistry_H
#define chemistry_H

#include "pointers.h"
#include <string>
#include <vector>
#include "mpi.h"
#include <fstream>
#include <sstream>
#include <math.h>       /* pow function */
/*
#include <stdlib.h>



#include <stdio.h>
#include <iostream>
#include <fstream>


#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
*/


#define MASTER 0

namespace MASKE_NS {
    
	class Chemistry : protected Pointers {
	public:
		
		Chemistry(class MASKE *);
		~Chemistry();
		
        std::vector<std::string> molnames;
        std::vector<double> mol_arad, mol_acir, mol_vapp, mol_ahyd, mol_bDH, mol_z;  // molecule linear sizes in radial and circumferential directions (when pertaining to a solid particle), apparent volume of molecule in solution, a and b parameters for Debye Huckel, molecule charge
	std::vector<int> mol_nufeb; // nufeb's chemical species index
	std::vector<int> mol_nufeb_form; // nufeb's chemical species form index: 0 = not hydrated, 1 = fully protonated, 2 = 1st deprotonated, 3 = 2nd deprotonated, 4 = 3rd deprotonated
	std::vector<std::string> mol_spec; // speciation element name  (I guess this is the name of the corresponding SOLUTION_MASTER_SPECIES in phreeqc)
        std::vector<double> mol_cins, mol_nins;  // vectors with concentration and number of molecules in solution: initially populated by solution.cpp
        std::vector<double> mol_cindV, mol_nindV;  // vectors with concentration and number of molecules in dV reservoire of solution (here assuming iso style of dV... vectors will be needed for aniso type)
        std::vector<std::string> gxnames, gxstyle;   // names and styles (const, or more complex ones possibly in future) of calculators of activity coeffs of activated complexes
        std::vector<std::string> DGnames, DGstyle;   // names and styles (const, or more complex ones possibly in future) of standard activation energy of activated complexes
        std::vector<double> cx, dim;   // vectors of standard state concentration and dimensionality of activated complexes associated to each DG
        std::vector<std::string> rxnames, chnames;  //name of reactions and chains
        std::vector<std::string> sennames, senstyle;  //name and style of surface energy calculators
        std::vector<std::string> chstyle;  //chain styles
        std::vector<std::string> mechnames, mechstyle, mechmode,mechinter;  //mechanism name, style, mode (straight or net), and energy scaling
        std::vector<int> rx_gxID,rx_DGID;   // pointers to gammax and DG calculators, specified for each simple reaction
        std::vector<double> Keq, ki;    // vectors of equilibrium constants and ki parameters (fraction contribution of interaction and surface energy) of chemical reactions
        std::vector<double> rx_dV_fgd, rx_dV_bkg;    // change of foreground and backgroun volume due to reaction
        std::vector<double> ch_dV_fgd, ch_dV_bkg;    // change of foreground and backgroun volume due to reaction chain
        std::vector<std::vector<double> > ch_rdV_fgd, ch_rdV_bkg;    // relative importance in terms of volme change for each step in chain reaction: wil be used in fix_delete sometimes to compute fraction of surface and volume change to be attributed to each reaction in chain
        std::vector<double> ch_arac;    // ratio between radial and circumferential size of chain
        std::vector<double> rx_ar_min, rx_ar_max, rx_ar_avg, rx_ar_avv, rx_ar_cum;    // change of foreground length in radial direction due to reaction: min, max, average (sum/N), volume-based average (sumV^1/3), cumulative (sum_ar)
        std::vector<double> rx_ac_min, rx_ac_max, rx_ac_avg, rx_ac_avv, rx_ac_cum;    // same as above, in circumferential direaction
        std::vector<bool> mechchain;   // vector saying whether each mechanisms refers to a chain of reactins (true) or to an individual reaction (false)
        std::vector<int> mechrcID, mechsenID;    // ID of reaction (or chain) and of surface energy calculator corresponding to each mechanism

        
        
        std::vector<std::vector<double> > gxcoef;    // a vector of vectors containing the coefficients for each gx calculator. How many coeffs are writtent and read depends on the style
        std::vector<std::vector<double> > DGcoef;    // a vector of vectors containing the coefficients for each DGx calculator. How many coeffs are writtent and read depends on the style
        std::vector<std::vector<double> > sencoef;    // a vector of vectors containing the coefficients for each surface energy calculator. How many coeffs are writtent and read depends on the style
        std::vector<std::vector<int> > bkg_molID,fgd_molID;    // a vector of vectors containing the IDs of molecules involved in each simple reaction.
        std::vector<std::vector<double> > bkg_nmol,fgd_nmol;    // a vector of vectors containing the stoichio coefficients for each simple reaction.
        std::vector<std::vector<bool> > bkg_isnotsolv;  // vector of vector contanining a flag saying whether each bkg molecule involved in reaction is or is not a solvent. If it is, the molecule will be excluded when computing supersaturations. true = is not a solvent, false = is a solvent
        std::vector<std::vector<int> > ch_rxID, ch_nrx;   // vectors of vectors containing the IDs and the number of reactions involved in each chain.

        
        

        
        int me;
        int Nmol, Nreax, Ngx,Nchain,NDG,Nsen,Nmech;
        
        void readDB(std::string);
        void addmolecule();
        void addgammax();
        void addDG();
        void addreax();
        void addsen();
        void addmech();
        void printall();
        double compDGx(int);
        double compgammax(int);
        double compSen(int);
		
	private:
        std::stringstream ss;

        
	};
	
}

#endif
