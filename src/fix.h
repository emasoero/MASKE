#ifndef FIX_H
#define FIX_H

#include "pointers.h"
#include <string>
#include <vector>
//#include "mpi.h"
//#include <string>
//#include <vector>
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#define MASTER 0

namespace MASKE_NS {
	
	class Fix : protected Pointers {
	public:
        Fix(class MASKE *);
		~Fix();
    
        int me;     // id of the current processor (rank)
        
        
        std::vector<std::string> allFixNames;   //a vector containing the name of all fixes irrespective of their type: used to check that each name is unique        
        
        std::vector<std::string> fKMCtype;   //type of rejection-free KMC event: can be delete, nucleate,..
        std::vector<std::string> fKMCname;   // user-defined name of KMC event
        std::vector<std::string> fKMCscom;   // subcomm in charge of attempting the event
        
        std::vector<std::string> fKMCmech;   // chemical mechanism to compute the rates of the event
        std::vector<std::string> fKMCreg;   // stored region name where to apply a fix nucleate
        std::vector<std::string> fKMClatt;   // stored lattice name to associate to a fix nucleate
        std::vector<std::string> fKMCmin;   // stored minimiser associated to a fix nucleate
        std::vector<int> fKMCmid;         // numeric id of mechanism in chemistry database
        std::vector<int> fKMCrid;         // numeric id of stored region
        std::vector<int> fKMClid;         // numeric id of stored lattice
        std::vector<int> fKMCminid;         // numeric id of stored minimiser

        std::vector<std::string> fKMCpgroup;    // name of lammps group to which a nucleated particle is added, besides "all"
        
        std::vector<int> fKMCptypeTRY;   // trial particle type used in fix nucleate
        std::vector<int> fKMCptype;  // particle type to which the event applies
        std::vector<std::string> fKMCpgeom;   // geometry type (e.g. sphere, ellipsoid) of trial particle to nucleate
        std::vector<double> fKMCpdiam;   // diameter of trial particle to nucleate
        
        std::vector<double> fKMC_DV;   // sampling volume associated with fix (useful only for nucleation fixes, and defined while sampling the events via a user-defined, stored command associated with the lattice, so here I just set it to zero
        
        std::vector<int> fKMCaID;         // ID of local fix in all-fix vector, which will enable any processor to retrieve fix-specific info later. Useful when a KMC event is accepted and to be executed
        std::vector<double> fKMCeveryt;         // time gap between two successive event trials
        std::vector<double> fKMCleval;         // time of last evaluation of process
        std::vector<int> fKMCnevents;         // number of events to sample for this process (dynamically updated in krun during the simulation)
        std::vector<int> fKMCglobID;         // position of local process in global list of all fKMC processes
        std::vector<int> fKMCrank;         // position of process in the local list of fKMC events of same type in this subcomm
        std::vector<int> fKMCfirst;         // position of first event in local rate etc vectors recoredr in fix_type.cpp
        std::vector<std::string> fKMCwei;   // type of weight calculator
        std::vector<std::vector<double> > fKMCwarg;    // arguments of weight calculator
        std::vector<std::string> fKMCsinST;   // stlye of input solution to compute beta (fixed or changing following the reactions?)
        std::vector<std::string> fKMCsinUL;   // input solution taken uniform or local
        std::vector<std::string> fKMCsinbox;   // input solution taken from box or box+dV
        std::vector<std::string> fKMCsoutUL;   // output solution taken uniform or local
        std::vector<std::string> fKMCsoutbox;   // output solution taken from box or box+dV


        std::vector<std::string> Ctype; // local-subcom types of continuous processes
        std::vector<std::string> Cname;  // local-subcom names of continuous processes
        std::vector<std::string> Cscom;     //subcom of each continuous process
        std::vector<double> Cdt;            // time increment of each cont proc
        std::vector<int> CglobID; //    position of local process in global list of all cont processes
        std::vector<double> Cleval; // vector with time of last exectution of this process



        
        std::vector<std::string> afKMCtype;   // same as above but gathering all events and used by MASTER only
        std::vector<std::string> afKMCname;
        std::vector<std::string> afKMCscom;
        
        std::vector<std::string> afKMCreg;
        std::vector<std::string> afKMClatt;
        std::vector<std::string> afKMCmin;
        std::vector<std::string> afKMCmech;
        std::vector<int> afKMCmid;
        std::vector<int> afKMCrid;
        std::vector<int> afKMClid;
        std::vector<int> afKMCminid;
        
        std::vector<double> afKMC_DV;
        
        std::vector<int> afKMCptypeTRY;
        std::vector<int> afKMCptype;
        std::vector<std::string> afKMCpgeom;
        std::vector<double> afKMCpdiam;
        
    
        std::vector<std::string> afKMCwei;
        std::vector<std::vector<double> > afKMCwarg;
        std::vector<std::string> afKMCsinST;
        std::vector<std::string> afKMCsinUL;
        std::vector<std::string> afKMCsinbox;
        std::vector<std::string> afKMCsoutUL;
        std::vector<std::string> afKMCsoutbox;
        std::vector<double> afKMCeveryt;
        std::vector<double> afKMCleval;
        std::vector<int> afKMCnevents;
        std::vector<double> afKMCcumRate;         // vector containing the cumulated rates of all events in each process separately
        
        std::vector<std::string> Contype;   // type of continuous processes: can be grow, dissolve, diffuse
        std::vector<std::string> Relxtype;  // type of particle relaxation: can be minimize, NVT, NPT ..
        
        
        std::vector<std::string> aCtype;
        std::vector<std::string> aCname;
        std::vector<std::string> aCscom;
        std::vector<double> aCdt;

        
        
        void add(std::string);
        void printall();
        
	private:
        void check_name_unique(std::string);

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
	};
	
}

#endif
