/*****************************************************************************
 * FILE:
 * DESCRIPTION:
	- first argument is the number of kinetic time steps
	- second argument is the number of nodes on which lammps should run (must be 1 until new developments...)
 *
 ****************************************************************************/
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
//#define  tr1::arraySIZE	16000000
//#define  tr1::arraySIZE	16
#define  MASTER		0
#define M_PI 3.14159265359
#define AN 6.02214129E+23   //   1/mol
#define kB 1.3806488E-23    //   J/K
#define R 8.3144617         //   J/(mol*K)
#define hP 6.626069E-34     //   J*s
#include <vector>
#include <fstream>
#include "iostream"
//#include <tr1/tr1::array>
#include <array>
#include <string>
#include <sstream>
#include <math.h>
#include <unistd.h>
#include <iomanip>      // std::setprecision
#include <sstream>      // std::stringstream, std::stringbuf
#include <sys/types.h>
#include <dirent.h>

#include <cassert>

#include "lammps.h"     // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
using namespace LAMMPS_NS;

#define PARALLEL

template <typename T>
std::string to_string(T value)
{
	std::ostringstream os ;
	os << value ;
	return os.str() ;
}


// DECLARE STRUCTURES (run in all processsors)

struct Atom2{
	// Atom2 types in the simulation
	std::string name;
	double mass;
};
struct MolSol{
	// molecule in solution: type, concentration, and properties.
	std::string name;
	std::vector<double> atinm;	//number of attoms of each species in the molecule
	double conc;	// number of molecules per unit volume
	double ninbox;  // at such concentration, how many such molecules are there in the simulation box?
	double Ei;	//internal energy due to Atom2 interacions, viz. lammps Denergy if the molecule is broken apart into single attoms
	double Es;	//solvation energy
	double q;	//charge.... AT SOME POINT WE WILL NEED TO CHECK ELECTRONEUTRALITY...
	//something about density of states...
};
int const MAX_n_shape_attr = 6;  // This specifies the maximum number of attributes that allow to describe the most complex particle shape in the simulation. E.g., 1 for shperical, 3 for ellipsoidal, etc..
struct PartType{
	// type of particles that can nucleate or exist initially (e.g. a substrate)
	std::string name;
	int ptid;	//NB: this tells something about the chemical nature of the particle. A more specific type for LAMMPS is then specified for each particle (see Particle structure below), in order to attribute the correct interaction. For example, if two particles of CSH type have different diameters, they can be considered for different interaction potentials in LAMMPS.
	std::vector<double> atins;	//number of attoms of each species in a molecule of solid
	double volmol;	/// molecular volume (actual only at nucleation. Further reactions can alter this value in single particles)
	double Eisp;	// internal energy per unit volume (in future, might be recomputed in every particle if they can be subjected to changes in chemistry, e.g. due to leaching or sequestration)
	double Essp;	// surface energy per unit surface (in future, it might be recomputed to account for an evolving solution chemistry, local environment, and possible chemical changes of single particles)
	double dens;	// density
	std::vector<int> RXnuc;	//list of id of possible chemical reactions for nucleation
	std::vector<int> RXgrow;	//list of id of possible chemical reactions for growth
	std::vector<int> RXdiss;	//list of id of possible chemical reactions for dissolution
	//lists for leaching and sequestration will be needed too. Solid-solid mass transfers are also left for the future..
	//something about density of states...
	std::vector<double> nuclsize;	//list of possible sizes of nuclei
	std::vector< std::array<double, 3> > nuclori;	//list of possible orientations of nuclei (here with three Euler angles, but could be converted to quaternions in future)
	std::vector<std::string> nuclshape; // list of possible shapes: names. 
	std::vector< std::array<double, MAX_n_shape_attr> > shape_attr; //list of shape attributes tr1::tr1::arrays
	std::vector<int> n_attr; // the number of significant attributes in any entry of the shape_attr vector. If entry 0 is a sphere, n_attr[0]=1, if entry 1 is an ellipsoid, n_attr[1]=3. This is usefule with MPI.
};
struct ChemRX{
	// chemical reaction: type, stoichiometry, and activation energy.
	// FOR NOW, THESE ARE ALL SOLUTION-SOLID RECTIONS. SOLID-SOLID and SOLUTION-SOLUTION REACTIONS WILL REQUIRE FURTHER IMPLEMENTATION
	std::string name;
	std::vector<double> mtosol;	//number of molecules of each species added to the solution
	std::vector<double> atosld;	//number of attoms of each species added to the solid
	double DEiso;	// activation energy of this reaction in isolated conditions
	double DVsld;	// change of volume of the solid induced by the chemical reaction
};
struct Particle{
	// the particles in the simulation
	int type;	//particle type id.. the same that will be used in lammps! Must be one of the ptid in parttype
	int typepos;	// this is computed as the actual position in the parttype vector, from 0 to its size-1
	std::vector<double> atinp;	//total number of attoms of each species in the particle
	double vol, surf, diam, mass, dens;
    int int_id; //id for the interactions in LAMMPS
	double x,y,z;
	double vx,vy,vz;
	double Ei;	// internal energy (in future, might be recomputed in every particle if they can be subjected to changes in chemistry, e.g. due to leaching or sequestration)
	double Esurf;	// surface energy (in future, it might be recomputed to account for an evolving solution chemistry, local environment, and possible chemical changes of single particles)
	double o1,o2,o3; //orientation angles. Might be converted to quaternions in future.
	std::string shape;  //e.g., Sphere, Ellipsoid, etc.. This will define the rules for the growth!
	int n_attr;		// number of relevant attributes for the shape description, i.e. relevant entries in the vector below
	double sh[MAX_n_shape_attr];	// shape descriptors
};
struct Interact{
    // the types of interactions between different particle types. These will be used to generate interaction id's
    // to be given to LAMMPS, or just to assign each particles a group id that is consistent with the LAMMPS input file and potential
    // type (or table)
    std::string name;   //e.g., Lennard-Jones , Mei , PolyMei , PairTable , Gay-Berne , etc..
    int nptypes;        // number of particle types involved in such interactions (2 for pairwise, 3 for three body, etc..)
    std::vector<int> ids;   // list of partycle type id
    int n_attr;         //the name will imply a certain number of attributes
    std::vector<double> attr; // the list of numerical values for the attributes
    /* For example:
        Lennard-Jones has two attributes
     */
};


int main (int argc, char *argv[])
{
    //check immediately that all arguments are passed:
    if (argc < 3){
        printf("ERROR: not enough input arguments from command line");
        exit(0);
    }
    clock_t start, end;
    start = clock();
    
    // DECLARE VARIABLES, tr1::arrayS, AND VECTORS (run in all processsors)
    int natsp;//atoi(argv[1]); ;  // number of Atom2ic species in solution+solids
	std::vector<Atom2> attoms;
	std::vector<MolSol> molsol;
	std::vector<ChemRX> chemRX;
	std::vector<PartType> parttype;
	std::vector<Particle> parts;
    int npart;//atoi(argv[2]); //number of particles
    double box[3];
	std::string c0fname = "conf0.init";
	std::vector< std::array<int, 5> > nucType;
    int NKtsteps =  atoi(argv[1]);  // total numnber of kinetic time steps to be performed
    int npL = atoi(argv[2]);        // number of processors that you want to use for each lammps instance (lower-level parallelization)
    //hereafter, stuff for the parallel implemenetation
	int   numtasks, taskid, rc , tag1=2, tag2=1, source;
	int *chunksize, *offset;
    MPI_Status status;
    //srand(12);
	
    double T = atof(argv[3]);
    double supersat = atof(argv[4]);
    
    // DECLARE FUNCTIONS (run in all procs)
	void read_chem(std::vector<Atom2> &attoms,std::vector<MolSol> &molsol, std::vector<PartType> &parttype, std::vector<ChemRX> &chemRX);
	void read_conf(std::string c0fname , double box[], std::vector<Particle> &parts);
	void read_Pchem(std::vector<Particle> &parts,std::vector<PartType> &parttype,std::vector<Atom2> &attoms);
	void print_input(std::vector<Atom2> &attoms,std::vector<MolSol> &molsol, std::vector<PartType> &parttype, std::vector<ChemRX> &chemRXdouble,double box[], std::vector<Particle> &parts);
    
    
	#ifdef PARALLEL
	// INIT MPI (run in all procs)
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);	//number of processors employed
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid); //gets the id of the current processor
    printf ("MPI task %d has started, numtasks= %d, NKtsteps= %d\n", taskid, numtasks, NKtsteps);
	#else
	taskid=0;
    numtasks=1;
	#endif
    
	bool error_input=false;	//error flag, to stop all processors if there is an error in the input
	std::vector<std::string> err_msg;	// all error messages recorded during the input phase, to be piped out at the end of it
	
    // READ INPUTS
	
	// IMPORT CHEMICAL INFO AND POSSIBLE NATURES OF NUCLEI (all processors record these)
	read_chem(attoms,molsol,parttype,chemRX);
	#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
	#endif
	


//if (taskid==MASTER){
//srand(12);}
//else{
//srand(taskid*11);
//}

	//check that each possible chemical reaction conserves mass
	if (taskid == MASTER){
		double count_attoms[attoms.size()];
		for (int i=0; i<chemRX.size(); i++) {
			for (int ii=0; ii<attoms.size(); ii++)
	 {count_attoms[ii]=0.;
//	  std::cout<<"ChemRX["<<i<<"]= "<<chemRX[i].name<< " attoms["<<ii<<"]= "<<attoms[ii].name<<std::endl;
	 }

			for (int j=0;j<molsol.size(); j++) {
				for (int k=0;k<attoms.size(); k++) {
					count_attoms[k] += chemRX[i].mtosol[j] * molsol[j].atinm[k];
//	  std::cout<<"chemRX= "<<chemRX[i].name<<" molsol["<<j<<"]= "<<molsol[j].name<< " attoms["<<k<<"]= "<<attoms[k].name<<" mtosol["<<j<<"]= "<<chemRX[i].mtosol[j]<<" count_attoms= "<<count_attoms[k]<<std::endl;
				}
			}
			for (int j=0;j<attoms.size(); j++) {
				count_attoms[j] += chemRX[i].atosld[j];
//	  std::cout<<"ChemRX["<<i<<"]= "<<chemRX[i].name<< " attoms["<<j<<"]= "<<attoms[j].name<<" atosld= "<<chemRX[i].atosld[j]<<std::endl;
			}
			for (int ii=0; ii<attoms.size(); ii++) {
				if(count_attoms[ii]!=0.){
					printf("\n ERROR: chemical reaction %s does not respect mass balance \n",chemRX[i].name.c_str() );
					err_msg.push_back("\n ERROR: chemical reaction " + chemRX[i].name + " does not respect mass balance \n");
					error_input = true;
				}
			}
		}
	}
	
	if (taskid == MASTER){
		
		// IMPORT THE FIRST LAMMPS-STYLE CONFIGURATION
        //read_conf(c0fname,box,parts);
        
		// assign to each particle the correct id to its type position in the particle types vector
		{
			for (int i=0; i<parts.size(); i++) {
				bool found = false;
				for (int j=0; j<parttype.size(); j++) {
					if (parts[i].type == parttype[j].ptid) {
						parts[i].typepos = j;
						found = true;
					}
				}
				if (!found) {
					std::string s1 = to_string(i);
					std::string s2 = to_string(parts[i].type);
					printf("\n ERROR: particle %i has type %i which does not match with any specified particle type.",i,parts[i].type);
					err_msg.push_back("\n ERROR: particle " + s1 + " has type " + s2 + "which does not match with any specified particle type.");
					error_input = true;
				}
			}
		}

		// IMPORT CHEMICAL INFO FOR THE INITIAL CONFIGURATION (extra info compared to what lammps handles)
		//read_Pchem(parts,parttype,attoms);
        
        // READ THE INTERACTION INPUT FILE
        
        // CONSTRUCT AN INTERACTION TABLE AMONG ALL POSSIBLE PAIRS (AND NTUPLES FOR MULTI-BODY INTERACTIONS)
        
        // GIVE TO EACH PARTICLE ITS CORRECT INTERACTIO ID
        // All this part on the interactions, I leave it for later... For now lets just use a small number of simple potentials..
		
		//print a whole heap of stuff
//		print_input(attoms,molsol,parttype,chemRX,box,parts);
    }
    
	
    // INPUT ERROR HANDLING
	#ifdef PARALLEL
	// THE ERROR FLAG IS SENT FROM THE MASTER TO THE SLAVES
	if (taskid == MASTER){
		for (int dest=1; dest<numtasks; dest++) {
			MPI_Send(&error_input, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
		}
	}
	if (taskid > MASTER) {
		source = MASTER;
		MPI_Recv(&error_input, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
	}
	#endif
	// QUIT IF THERE IS AN ERROR
	if (error_input) {
		if (taskid == MASTER){
			for (int i=0; i<err_msg.size(); i++) {
				printf("%s",err_msg[i].c_str());
			}
		}
		#ifdef PARALLEL
		MPI_Finalize();
		#endif
		exit(0);
	}
    
    // THE MASTER PREPARES THE FIRST CONFIGURATION INPUT FILE FOR LAMMPS (e.g. data.init)
	if (taskid == MASTER){
        for (int i=0; i<parts.size(); i++) {
            //record each particle's position parts[i].x parts[i].y, velocity parts[i].vx parts[i].vy, and whatever needed, to the input config file for LAMMPS
            // ...to be written... for now I just use a data.init file that IGOR gave me and that I have put already into the RUN folder
        }
    }
    
    // IF NEEDED, THE MASTER PREPARES THE INTERACTION TABLE FILE FOR LAMMPS (unless someone has already prepared it consistently before, our of this code)
    if (taskid == MASTER){
    }
    
	// INIT HIGHER-LEVEL PARALLELIZATION (each processor creates its own working space)
    
    // npH SLAVE processors are used to start LAMMPS and extract+elaborate its outputs when possible
    // Here I tell each processor to create its own folder to run LAMMPS and to get a copy of the initial input configuration
    //	and input file for LAMMPS
    //
    //  NB: There is also a lower level parallelization which comes from the fact that each processor can run LAMMPS in parallel
    //		on npL cores, where npL stands for number of "SUPERSLAVE" processors used by each SLAVE processor. Overall, the number
    //		of cores used in a simulation in npH * npL
    //
    // Handling this 2-levels parallelization should not be difficult, but some issues may rise in clusters with queues.
    
    // Each Hlevel processor creates its working folder, and copies the LAMMPS executable and initial particle config there there
    std::string stskid = to_string(taskid);
    std::string workFold;		workFold = "tempLAMMPS_"+stskid;
    
    //system(("mkdir "+workFold).c_str());
    //system(("cp lmp_* "+workFold).c_str());
    //system(("cp data.init "+workFold).c_str());
    
    // if needed, copy the interaction table file too..
    std::string snpL = to_string(npL);  //converts number of lower-level processors into a string, to then use to call LAMMPS
	
	//INIT NUCLEATION (the master shares the job load)

	/*
	 // This records a list of id's that are used in the nucleation
	 // The id's refer to:
	 // - parttype
	 // - RXnuc in parttype
	 // - nuclsize in parttype
	 // - nuclshape in parttype
	 // The result is a vector of tr1::array[4]. The length of the vector is N, with:
	 //		N = SUM_i  [    parttype[i].RXnuc.size() * parttype[i].nuclsize.size() * parttype[i].nuclshape.size()   ]
	 // i.e.   (number of possible chemical reaction)   *   (number of possible nuclei sizes)  *    (number of possible nuclei shapes) 
	 //			summed over all the possible particle types
	 //
	 // N is divided in npH chunks, where npH is the number of processors involve in the high-level parallelization (see above)
	 // Each chunk is passed to a processor
	 //
	 */
		
	int Nacc[numtasks];     //an tr1::array of number of nucleations accepted during the generic iteration
	int acc_size=0.;		// total number of accepted nucleations during the generic iteration
	
	// record the id's of parttype, RXnuc in parttype, nuclsize in parttype, and nuclshape in parttype
	for (int i=0; i<parttype.size(); i++) {
	  for (int j=0; j<parttype[i].RXnuc.size(); j++) {
	    for (int k=0; k<parttype[i].nuclsize.size(); k++) {
              for (int l=0; l<parttype[i].nuclori.size(); l++) {
                for (int m=0; m<parttype[i].nuclshape.size(); m++) {
                        std::array<int, 5> temp_type;
                        temp_type[0] = i;
                        temp_type[1] = j;
                        temp_type[2] = k;
                        temp_type[3] = l;
                        temp_type[4] = m;
                        nucType.push_back(temp_type);
//std::cout<<"nucType: ptype.size_i="<<i<<" RXnuc_j="<<j<<" size_k="<<k<<" ori_l="<<l<<" shape_m="<<m<<" nucType.size="<<nucType.size()<<std::endl;
                    }
                }
            }
        }
    }
	//subdivide the nucType vector in npH chunks, such that each node gets as few jobs as possible
	chunksize = new int[numtasks];
	for (int i=0; i<numtasks; i++) {chunksize[i]=0;}
	
	{
		bool finito = false;
		chunksize[0] = nucType.size();	//first assign all nucleation jobs to the master
		while (!finito) {
			finito = true;
			// for all the slaves (i starts from 1!)
			for (int i=1; i<numtasks; i++) {
				if (chunksize[0]-chunksize[i] > 1) {
					chunksize[0]--;
					chunksize[i]++;
					finito = false;
                }
			}
		}
//for (int i=0; i<numtasks; i++) {std::cout<<"After equilibration chunksize["<<i<<"]= "<<chunksize[i]<<std::endl;}
	}

	// set offsets for each processor to read its own nucleation jobs
	offset = new int[numtasks];
	for (int i=0; i<numtasks; i++) offset[i]=0;
	
	//master's offset is 0, while for the slaves (i starts from 1) is...
	for (int i=1; i<numtasks; i++) offset[i] = offset[i-1] + chunksize[i-1];

for(int i=0; i<numtasks; i++)
std::cout<<" offset["<<i<<"]= "<<offset[i]<<" chunksize["<<i<<"]= "<<chunksize[i]<<std::endl;	
	
	// All processors use their part of the nucType vector to create their lammps input scripts for nucleation.
	// For now, I just copy the lj.in to each working folder and assign a tentative lattice size of 100
	//system(("cp in.lj "+ workFold).c_str());
	int ntry = 100;		//this should be either given consistently with LAMMPS' lattice, or read from LAMMPS output file
	//create a number of vectors that will be used later in the nucleation routine

   	std::vector<double> Utry, x_try, y_try, z_try, diam_try, o1_try, o2_try, o3_try;
	std::vector<int> ptype_try, shape_try;
	std::vector<int> ptype_acc, shape_acc;
	std::vector<double> x_acc, y_acc, z_acc, diam_acc, o1_acc, o2_acc, o3_acc;
	std::vector<int>acc; std::vector<int>acc_chunk;   // added IS	
 
    // START MAIN LOOP
    // F&F: initiate lammps here, which will be kept running for the whole chunk
    // BEWARE!!!!! WE ARE OPENING LAMMPS IN ALL PROCESSORS, NOT ONLY ON HIGH-LEVEL PARALLELISATION ONES. THIS ASSUMES IMPLICITLY THAT THERE IS NO LOWER LEVEL PARALLELISATION, IE LAMMPS RUNS IN SERIES!

    ptype_acc.clear(); shape_acc.clear(); x_acc.clear(); y_acc.clear(); z_acc.clear(); diam_acc.clear(); o1_acc.clear(); o2_acc.clear(); o3_acc.clear(); acc.clear(); acc_chunk.clear();
    int count_acc = 0;
    
	MPI_Barrier(MPI_COMM_WORLD);
    
        LAMMPS *lmp;
        lmp = new LAMMPS(0,NULL,MPI_COMM_WORLD);
        double R_insert = 0.0 , R_delete = 0.0;         // individual insert & delete rates
        int Ktstep = 0, inserted = 0, deleted = 0;
        int prev_Ktstep = 0, prev_inserted = 0, prev_deleted = 0, prev_N2 = 0;
        double dt_discard = 0., KMCtime = 0., prev_ins = 0., prev_del = 0., prev_tot = 0.;
    
        //double gamma = 8.8e-20;                         // Nonat homogeneous [J/nm2]        6.05e-21  - hetero
        double nu = 1.;                                 // rate exp prefactor
        //int n = ceil((4.*M_PI*pow(2.5,3)/(3.*0.26768))-0.5);                // 0.26768 nm3   volume of CSHII molecule Bullard
        //double omega = 4.*M_PI*pow(2.5,2);
        //double alpha_c = -log(supersat)+gamma*omega/(n*kB*T);
        //double alpha_m = 1/(2*n*kB*AN*T);
        std::vector<double> time_list;
        std::vector<int> Nucl_list;
        std::vector<std::string> content;
        std::ofstream outhistory, outdump, outtest;
        std::ifstream ifhistory("history.dat"), ifdump("BOX.dump");
    
    if (taskid == MASTER){                              /* Setting up restart positions */
        outtest.open("test.txt");
        std::string line;
        std::stringstream lastline;
        
        if (!ifhistory) {
            outhistory.open("history.dat");
            outhistory<<"# KMCtstep \t +N2 \t -N2 \t N2total \t Timestep \t Cumul_Time \t Cumul_R+ \t Cumul_R- \t Cumul_R"<<std::endl;
            outhistory<<"# T="<<T<<"\t beta="<<supersat<<std::endl;
            Ktstep = 0, prev_inserted = 0, prev_deleted = 0, prev_N2 =0;
            dt_discard = 0.0, KMCtime = 0.0, prev_ins = 0.0, prev_del = 0.0, prev_tot = 0.0;
        }
        else {
            while (std::getline(ifhistory, line)) content.push_back(line);
            
            outhistory.open("history.dat", std::ios::app);
            lastline << content[content.size()-1];
            lastline >> prev_Ktstep >> prev_inserted >> prev_deleted >> prev_N2 >> dt_discard >> KMCtime >> prev_ins >> prev_del >> prev_tot;

            prev_Ktstep = prev_Ktstep + 1;
            prev_inserted = prev_inserted + 1;
            if (prev_N2 == 0) prev_N2 = 0;
            else prev_N2 = 1;
        }
        
        if (!ifdump) {outdump.open("BOX.dump");}
        else {outdump.open("BOX.dump", std::ios::app);}
    }

    DIR *dp;                                        /*  Choosing what data.init file to read */
    struct dirent *ep;
    dp = opendir ("./Data/");
    char src[200] , max[200]="data.init" ;
    
    if (dp != NULL)
    {
        while ((ep = readdir (dp))){
            strcpy( src , ep->d_name);
            //puts (src);
            if (strcmp(max , src) < 0) strcpy(max,src);
        }
    }
    (void) closedir (dp);
    
        time_t t;
        srand((unsigned) time(&t));
        rand();
        
        lmp->input->one("units           real");
        lmp->input->one("atom_style      sphere");
        lmp->input->one("boundary        f f f");
        lmp->input->one("pair_style      lj/cut 110");
    
        std::string data_init = max;
        std::stringstream ss_data;
        ss_data<<"read_data ./Data/"<<data_init;
        //lmp->input->one((ss_data.str()).c_str());
        lmp->input->one("timestep        1");
        lmp->input->one("region          Sim_Box block 0 1600 0 1600 0 1600");
        lmp->input->one("create_box      4 Sim_Box");
    
        lmp->input->one("pair_coeff      * * 0 50");
        lmp->input->one("pair_coeff      1 1 0 50");
        lmp->input->one("pair_coeff      2 2 500 50");
        lmp->input->one("pair_coeff      3 3 0 50");
        lmp->input->one("pair_coeff      1 2 500 50");
        lmp->input->one("pair_coeff      1 3 500 50");
        lmp->input->one("pair_coeff      2 3 500 50");
    
        lmp->input->one("lattice         sc 50 origin 0 0 0 orient x  1 0 0 orient y  0 1 0 orient z  0 0 1");
        lmp->input->one("region          Layer sphere 800 800 800 25 units box");
        lmp->input->one("create_atoms    1 region Layer");
        lmp->input->one("lattice         sc 50 origin 0 0 0 orient x  1 0 0 orient y  0 1 0 orient z  0 0 1");
        lmp->input->one("region          SphBox sphere 800 800 800 25 side out units box");
        lmp->input->one("region          OutBox block 100 1500 100 1500 100 1500 units box");
        lmp->input->one("region          Box intersect 2 OutBox SphBox");
        lmp->input->one("create_atoms    3 region Box");
    
        lmp->input->one("group           Layer type 1");
        //lmp->input->one("group           Nuclei type 2");
        lmp->input->one("group           Trial_Particles type 3");
        lmp->input->one("set             group Layer diameter 50");
        lmp->input->one("set             group Trial_Particles diameter 50");
        lmp->input->one("variable        n50 atom 100");
        lmp->input->one("set             group all mass v_n50");
    
        //lmp->input->one("fix             Walls all wall/reflect zlo EDGE zhi EDGE");
        lmp->input->one("fix             Freeze_1 Layer setforce 0.0 0.0 0.0");
    
        lmp->input->one("neighbor        2.0 bin");
        lmp->input->one("neigh_modify    exclude group Layer Layer");
        lmp->input->one("neigh_modify    exclude group Trial_Particles Trial_Particles");
        lmp->input->one("neigh_modify    every 2 delay 0 check yes page 100000 one 10000");
    
        lmp->input->one("compute         PE all pe/atom");
        lmp->input->one("compute         Radius all property/atom radius");
        lmp->input->one("compute         Temp all temp");
        lmp->input->one("compute_modify  Temp dynamic yes");
        lmp->input->one("variable        N1 equal count(Layer)");
        lmp->input->one("variable        N2 equal count(Nuclei)");
        lmp->input->one("variable        N3 equal count(Trial_Particles)");
    
        lmp->input->one("thermo_style    custom step temp pe etotal press vol v_N2");
        lmp->input->one("thermo_modify   lost warn flush yes");
        lmp->input->one("thermo_modify   temp Temp");
        lmp->input->one("thermo          200");
        //lmp->input->one("dump            PE all custom 100 PE.dump id type diameter c_PE x y z");
        //lmp->input->one("dump_modify     PE sort id");
    
        for (Ktstep; Ktstep < NKtsteps; Ktstep++) {
        if (taskid == MASTER){
            std::cout<<"           ##############################################################"<<std::endl;
            std::cout<<"           ####                      Ktstep = "<<prev_Ktstep + Ktstep<<"                      ####"<<std::endl;
            std::cout<<"           ##############################################################"<<std::endl;
        }
            
            lmp->input->one("group           Nuclei type 2");                   //indentified here so v_N2 can count new nuclei
            lmp->input->one("group           Trial_Particles type 3");
            lmp->input->one("set             group Nuclei diameter 50");
            
            std::stringstream ss1;
            ss1<<std::setprecision(10)<<std::fixed<<"displace_atoms  Trial_Particles random 1 1 1 "<<rand()<<" units box";
            lmp->input->one((ss1.str()).c_str());
            
            lmp->input->one("fix             Freeze_2 Nuclei setforce 0.0 0.0 0.0");
            lmp->input->one("pair_coeff      2 2 0 50");
            lmp->input->one("pair_coeff      1 2 0 50");
            lmp->input->one("min_style       quickmin");
            lmp->input->one("min_modify      dmax 0.1 25 25 25");
            lmp->input->one("minimize        0.0 0.0 1000 1000");
            //lmp->input->one("delete_atoms    overlap 40 Trial_Particles Nuclei");
            lmp->input->one("set             group all mass v_n50");
            lmp->input->one("run             0");
            
            // unsorted energies and ids of all particles
            int nlocal = static_cast<int> (lmp->atom->nlocal);
            int *aID = new int[nlocal];
            aID = ((int *) lammps_extract_atom(lmp,"id"));
            double *E = new double[nlocal];
            E = ((double *) lammps_extract_compute(lmp,"PE",1,1));     // local E
            int *aTYPE = new int[nlocal];
            aTYPE = ((int *) lammps_extract_atom(lmp,"type"));
            double *Rad = new double[nlocal];
            Rad = ((double *) lammps_extract_compute(lmp,"Radius",1,1));
            
//##################### Assembling unsorted aID vectors on MASTER ###########################
            int natoms = static_cast<int> (lmp->atom->natoms);
            int aIDglob[natoms];
            double Eglob[natoms];
            double Rglob[natoms];
            int nlocs[numtasks];
            
            nlocs[taskid] = nlocal;

            if (taskid > MASTER){
                int dest = MASTER;
                MPI_Send(&nlocs[taskid], 1, MPI_INT, dest, tag2, MPI_COMM_WORLD);
            }
            if (taskid == MASTER) {
                for (int i=1; i<numtasks; i++) {
                    int source = i;
                    MPI_Recv(&nlocs[i], 1, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                }
            }
            MPI_Bcast(&nlocs, numtasks, MPI_INT, MASTER, MPI_COMM_WORLD);

            // find the right position where the current processor should record its local ids in the global id vector
            int right_pos=0;
            for (int ii=0; ii<taskid; ii++) {right_pos += nlocs[ii];}
           
            // record local ids in global id vecotr (before communicating it to other processors)
            for (int ii=0; ii<nlocal; ii++) {
                aIDglob[right_pos+ii] = aID[ii];
                Eglob[right_pos+ii] = E[ii];
                Rglob[right_pos+ii] = Rad[ii];
            }
            
            int tag3 = 3;
            // send global id vectors to the master
            if (taskid > MASTER){
                int dest = MASTER;
                MPI_Send(&aIDglob[right_pos], nlocal, MPI_INT, dest, tag2, MPI_COMM_WORLD);
                MPI_Send(&Eglob[right_pos], nlocal, MPI_DOUBLE, dest, tag1, MPI_COMM_WORLD);
                MPI_Send(&Rglob[right_pos], nlocal, MPI_DOUBLE, dest, tag3, MPI_COMM_WORLD);
            }
            
            if (taskid == MASTER) {
                for (int i=1; i<numtasks; i++) {
                    right_pos = 0;
                    for (int ii=0; ii<i; ii++) {
                        right_pos += nlocs[ii];
                    }
                    int source = i;
                    MPI_Recv(&aIDglob[right_pos], nlocs[i], MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&Eglob[right_pos], nlocs[i], MPI_DOUBLE, source, tag1, MPI_COMM_WORLD, &status);
                    MPI_Recv(&Rglob[right_pos], nlocs[i], MPI_DOUBLE, source, tag3, MPI_COMM_WORLD, &status);
                }
            }   // The master now should know all the unusorted ids aIDglob
            MPI_Barrier(MPI_COMM_WORLD);
//##################### End of Assembly #####################################################

            //sorted positions, ids, and types
            double *x = new double[3*natoms];
            lammps_gather_atoms(lmp,"x",1,3,x);
            int *aIDs= new int[natoms];
            lammps_gather_atoms(lmp,"id",0,1,aIDs);
            int *aTYPEs = new int[natoms];
            lammps_gather_atoms(lmp,"type",0,1,aTYPEs);
            
            // generate a vector of ids to go from sorted to unsorted vectors
            int *id2id = new int[natoms];     //per-atom E is in unsorted[natoms] while id & types are in sorted[3*natoms] so we must match them
            
            if (taskid == MASTER) {
                for (int ii=0; ii<natoms; ii++) {
                    for (int kk=0; kk<natoms; kk++) {
                        if (aIDglob[kk]==aIDs[ii]) {
                            id2id[ii]=kk;
                        }
                    }
                }
            }
            
            double tN1 = *((double *) lammps_extract_variable(lmp,"N1",0));
            int N1 = (int) tN1;
            double tN2 = *((double *) lammps_extract_variable(lmp,"N2",0));
            int N2 = (int) tN2;
            double tN3 = *((double *) lammps_extract_variable(lmp,"N3",0));
            int N3 = (int) tN3;
             
            //list of rates for all trial particles and deleted nuclei
            std::vector<int> pins;      //pointing to the position of current trial particle in sorted vectors
            std::vector<double> rates;  //rates for trial particle only (consist of rate_insert and rate_delete elements)
            std::vector<double> rate_insert;  //insertion rates for trial particle only
            std::vector<double> rate_delete;  //deletion rates for trial particle only
            std::vector<double> Mtype;  //type of move associated with the rate (0 = insertion, 1 = deletion, 2 = growth)
            
            //######## INSERTION ####### WE ARE FILLING RATES VECTOR FOR EACH TRIAL PARTICLE (POSSIBLE INSERTION 3 -> 2)
            double Vc = pow(50.,3);                         // V of Cell
            double avN2 = N2/(1E+9/Vc);
            double phi = avN2*4.*M_PI*pow(25.,3)/(3.*Vc);
            double a3 = pow(50.,3);                         // V[A^3] of CSHII molecule

            if (taskid==MASTER) {
                //for all particles (in sorted vector)
                for (int ii=0; ii<natoms; ii++) {
                    if (aTYPEs[ii]==3) {
                        pins.push_back(ii);
                        Mtype.push_back(0);
                        R_insert = nu*(Vc/a3)*exp(-2*E[id2id[ii]]/T + supersat/*-log((1.-phi)*Vc/a3)*/);
                        rates.push_back(R_insert);
                        rate_insert.push_back(R_insert);
                    }
                    else if (aTYPEs[ii]==2) rate_insert.push_back(0);
                }
            }

            //######## DELETION ####### CYCLE OVER ALL THE EXISTING PARTICLES OF TYPE 2 AND TRY TO DELETE THEM TO ADD THE ASSOCIATED DELETION RATES
        if (taskid == MASTER) std::cout<<"!!!!!!!!!!!!!! DELETABLE NUCLEI = "<<N2<<std::endl;
            
            lmp->input->one("pair_coeff   2 3 0 50");
            lmp->input->one("pair_coeff   1 2 500 50");
            lmp->input->one("pair_coeff   2 2 500 50");
            lmp->input->one("run	  0");
            double *Enew = new double[nlocal];
            Enew = ((double *) lammps_extract_compute(lmp,"PE",1,1));     // local Enew
            int *aIDnew = new int[nlocal];
            aIDnew = ((int *) lammps_extract_atom(lmp,"id"));

//##################### Assembling NEW unsorted Enewglob vectors on MASTER ###########################
            int aIDnewglob[natoms];
            double Enewglob[natoms];
            
            nlocs[taskid] = nlocal;
  
            // find the right position where the current processor should record its local ids in the global id vector
            right_pos=0;
            for (int ii=0; ii<taskid; ii++) {right_pos += nlocs[ii];}
            
            for (int ii=0; ii<nlocal; ii++) {
                aIDnewglob[right_pos+ii] = aIDnew[ii];
                Enewglob[right_pos+ii] = Enew[ii];
            }

            if (taskid > MASTER){
                int dest = MASTER;
                MPI_Send(&aIDnewglob[right_pos], nlocal, MPI_INT, dest, tag2, MPI_COMM_WORLD);
                MPI_Send(&Enewglob[right_pos], nlocal, MPI_DOUBLE, dest, tag1, MPI_COMM_WORLD);
            }

            if (taskid == MASTER) {
                for (int i=1; i<numtasks; i++) {
                    right_pos = 0;
                    for (int ii=0; ii<i; ii++) {
                        right_pos += nlocs[ii];
                    }
                    int source = i;
                    MPI_Recv(&aIDnewglob[right_pos], nlocs[i], MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&Enewglob[right_pos], nlocs[i], MPI_DOUBLE, source, tag1, MPI_COMM_WORLD, &status);
                }
            }   // The master now should know NEW usorted Enewglob
            MPI_Barrier(MPI_COMM_WORLD);
//##################### End of Assembly #####################################################
            
            int PPos = 0, pos;
            double dt;
            if (taskid == MASTER){
                for (int ii=0; ii<natoms; ii++) {
                    if (aIDnewglob[ii] != aIDglob[ii]) {
                        std::cout<<"ERROR WITH THE IDS. MORE WORK NEEDED"<<std::endl;
                        exit(0);
                    }
                }

                for (int ii=0; ii<natoms; ii++) {
                    if (aTYPEs[ii]==2) {
                        pins.push_back(ii);
                        Mtype.push_back(1);
                        R_delete = nu*exp(2*Enewglob[id2id[ii]]/T);
                        rates.push_back(R_delete);
                        rate_delete.push_back(R_delete);
                    }
                    else if (aTYPEs[ii]==3) rate_delete.push_back(0);
                }
                //######### END of DELETION #################
                
                //turn rate vector into a CUMULATIVE one
                for (int ii=1; ii<rates.size(); ii++) {
                    rates[ii]+=rates[ii-1];
                    rate_insert[ii]+=rate_insert[ii-1];
                    rate_delete[ii]+=rate_delete[ii-1];
                    // if(ii >= 998) std::cout<<ii<<" atomic_id="<<aIDs[pins[ii]]<<"  ratE="<<std::setprecision(10)<<rates[ii]<<"  E="<<std::setprecision(10)<<E[pins[ii]]<<"  Mtype="<<Mtype[ii]<<" x="<<x[3*pins[ii]]<<" y="<<x[3*pins[ii]+1]<<" z="<<x[3*pins[ii]+2]<<std::endl;
                }
                
                //######### BINARY SEARCH ################# find event with chosen cumulative Rate
                double chosenR =  (double)rand() / (double)RAND_MAX *rates[rates.size()-1];
                int pre = 0;
                int post=rates.size()-1;
                
                while (pre < post) {
                    pos = (int)((pre+post)/2.);
                    if (rates[pos] < chosenR) pre=pos+1;
                    else post=pos;
                }
                pos=pre;
                
                std::cout<<" ########## ATTEMPTED NUCLEI Ktstep= "<<Ktstep<<" pos="<<pos<<"  atomic_id="<<aIDs[pins[pos]]<<" MType="<<Mtype[pos]<<" ratE="<<std::setprecision(10)<<rates[pos]-rates[pos-1]<<" x="<<x[3*pins[pos]]<<" y="<<x[3*pins[pos]+1]<<" z="<<x[3*pins[pos]+2]<<std::endl;
                //######### END OF SEARCH #################
                
                //lmp->input->one("group           Nuclei type 2");
                
                /* ######################## WRITING BOX and HISTORY files ############################## */
                    Nucl_list.push_back(N2);
                    time_list.push_back(KMCtime);
                    
                    dt = 1./rates[rates.size()-1]*log((double)RAND_MAX / (double)rand());
                    KMCtime += dt;
                    double rad2 = pow((x[3*pins[pos]] - 800),2) + pow((x[3*pins[pos]+1] - 800),2) + pow((x[3*pins[pos]+2] - 800),2);
                
                    outhistory<<prev_Ktstep + Ktstep<<"\t"<<prev_inserted + inserted<<"\t"<<prev_deleted + deleted<<"\t"<<prev_N2 + N2<<"\t"<<std::setprecision(40)<<dt<<"\t"<<KMCtime<<"\t"<<std::setprecision(10)<<prev_ins + rate_insert[rate_insert.size()-1]<<"\t"<<prev_del + rate_delete[rate_delete.size()-1]<<"\t"<<prev_tot + rates[rates.size()-1]<<"\t"<<sqrt(rad2)/*supersat+log((1.-phi)*Vc/a3)*/<<std::endl;
                    
                    outdump<<N1 + N2<<std::endl;
                    outdump<<"id   type   radius    x    y    z"<<std::endl;
                    for (int i = 0; i < natoms; i++)
                        if (aTYPEs[i] < 3) outdump<<aIDs[i]<<" "<<aTYPEs[i]<<" "<<Rglob[id2id[i]]<<" "<<x[3*i]<<" "<<x[3*i+1]<<" "<<x[3*i+2]<<std::endl;
                
                PPos = pins[pos];
            }
            
            if (taskid == MASTER){          // Send-Receive positions in x-array for chosen particle
                for (int dest=1; dest<numtasks; dest++) {
                    MPI_Send(&PPos, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
                }
            }
            if (taskid > MASTER) {
                source = MASTER;
                MPI_Recv(&PPos, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            }
            
            int Mtype_sr;                   // Mtype send-receive
            
            if (taskid == MASTER) {
                Mtype_sr = Mtype[pos];

                for (int dest=1; dest<numtasks; dest++) {
                    MPI_Send(&Mtype_sr, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
                }
            }
            if (taskid > MASTER) {          // Send-Receive chosen particle M-type
                source = MASTER;
                MPI_Recv(&Mtype_sr, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            }
            
            if (Mtype_sr==0) {
                //ACTUALL NUCLEI INSTERION on every processor. Due to domain decomposit apparently only one processor inserts/delets
                std::stringstream ss;
                ss<<std::setprecision(10)<<std::fixed<<"create_atoms 2 single "<<x[3*PPos]<<" "<<x[3*PPos+1]<<" "<<x[3*PPos+2]<<" units box";
                if (taskid == MASTER) std::cout<<" ########## ATTENTION INSERTION CHOSEN "<<ss.str()<<std::endl;

                lmp->input->one((ss.str()).c_str());
                lmp->input->one("set          group Nuclei diameter 50");
                inserted += 1;
            }
            
            else if(Mtype_sr==1){
                //ACTUAL NUCLEI DELETION on every processor
                std::stringstream ss;
                ss<<"set atom "<<aIDs[PPos]<<" type 4";
                if (taskid == MASTER) std::cout<<" ########## ATTENTION DELETION CHOSEN pos="<<pos<<" atomic_id="<<aIDs[PPos]<<" "<<ss.str()<<std::endl;
            
                lmp->input->one((ss.str()).c_str());
                lmp->input->one("group        Deletable type 4");
                lmp->input->one("delete_atoms group Deletable");
                deleted += 1;
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            lmp->input->one("delete_atoms group Trial_Particles");
            lmp->input->one("unfix            Freeze_2");
            lmp->input->one("set              group all mass v_n50");
            lmp->input->one("min_style        cg");
            lmp->input->one("min_modify       dmax 0.1 25 25 25");
            lmp->input->one("minimize         0.0 0.0 250 250");
            
            std::stringstream ss_restart;
            ss_restart<<"write_data       ./Data/data.init."<<prev_Ktstep + Ktstep<<".*";
            if((prev_Ktstep + Ktstep) % 50 == 0) lmp->input->one((ss_restart.str()).c_str());
            
            lmp->input->one("create_atoms     3 region Box");
            lmp->input->one("set              group all mass v_n50");
            lmp->input->one("pair_coeff       1 2 500 50");
            lmp->input->one("pair_coeff       1 3 500 50");
            lmp->input->one("pair_coeff       2 3 500 50");
            lmp->input->one("pair_coeff       2 3 500 50");
            lmp->input->one("pair_coeff       3 3 0 50");
            lmp->input->one("run              0");
        }   /* End of Ktstep loop */
    
        delete lmp;
        end = clock();
    
    if (taskid == MASTER){
        outhistory <<"Execution time: "<< (double)(end-start)/CLOCKS_PER_SEC <<" seconds "<<(double)(end-start)*(1/3600.)/CLOCKS_PER_SEC<<" hours"<< std::endl;
        outhistory.close();
        outdump.close();
        outtest.close();
        printf("\nSUCCESS\n\n");
    }
    
    #ifdef PARALLEL
    MPI_Finalize();
    #endif

    

/*
            // KINETIC CRITERION ONLY ON ATOMS OF TYPE 3
        
             Utry.clear(); x_try.clear(); y_try.clear(); z_try.clear(); diam_try.clear(); o1_try.clear(); o2_try.clear(); o3_try.clear(); ptype_try.clear(); shape_try.clear();
			
			// each processor reads the last entry of the LAMMPS output energy and xyz files, and records them to temporary vectors
			//  Igor, you should already have this
			for (int i=0; i<ntry; i++) {

           		Utry.push_back( - (double)rand()/(double)RAND_MAX);
				x_try.push_back( (double)rand()/(double)RAND_MAX);
				y_try.push_back( (double)rand()/(double)RAND_MAX);
				z_try.push_back( (double)rand()/(double)RAND_MAX);
				diam_try.push_back( (double)rand()/(double)RAND_MAX);	// this can be taken from nucType vector
				ptype_try.push_back(0);   // all CSH					// this can be taken from nucType vector
				shape_try.push_back(0);   // all spheres				// this can be taken from nucType vector
				o1_try.push_back( (double)rand()/(double)RAND_MAX);
				o2_try.push_back( (double)rand()/(double)RAND_MAX);
				o3_try.push_back( (double)rand()/(double)RAND_MAX);

//std::cout<<"PROC "<<taskid<<" chunk["<<nck<<"]"<<" Utry=["<<i<<"]= "<<Utry[i]<<" ptype= "<<ptype_try[i]<<" x= "<<x_try[i]<<" diam="<<diam_try[i]<<std::endl;
			}
			
			// each processor chooses which new trial nuclei satisfy the criterion to become real.
			// For each accepted nucleus, several vectors are created, including:
			//  xpos    ypos    zpos     diameter     particle_type     shape_id     orientation1   o2   o3
			//std::vector<int> id_acc;
			for (int i=0; i<ntry; i++) {
				if (Utry[i] < -0.99) {   //HERE WE SHOULD STICK THE KINETIC CRITERION!!    That's where the info in nucType become useful.
							// For now, I simply accept all particles with energy less than a certain value
			//id_acc.push_back(i);

            x_acc.push_back(x_try[i]);
			y_acc.push_back(y_try[i]);
			z_acc.push_back(z_try[i]);
			diam_acc.push_back(diam_try[i]);
			ptype_acc.push_back(ptype_try[i]);
			shape_acc.push_back(shape_try[i]);
			o1_acc.push_back(o1_try[i]);
			o2_acc.push_back(o2_try[i]);
			o3_acc.push_back(o3_try[i]);
            count_acc++;
                    
std::cout<<"Acc loop PROC["<<taskid<<"] chunk "<<nck<<" Utry["<<i<<"]= "<<Utry[i]<<" ptype= "<<ptype_acc[count_acc - 1]<<" x="<<x_acc[count_acc - 1]<<" diam= "<< diam_acc[count_acc - 1]<<" Acc in chunk["<<nck<<"] = "<<count_acc<<std::endl;
				}
			}
          
            acc.push_back(count_acc);
//std::cout<<"count acc= "<<acc[nck]<<std::endl;
           
           //std::cout<<" Proc "<<taskid<<" is HERE"<<std::endl;
		}  // end of lammps-create-accept loop
        
        #ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);	 //wait for all processors to have finished with their part
        #endif
        
        //each processor populates its part of the tr1::array containing the number accepted nucleations
        acc_size = count_acc ;
        Nacc[taskid] = acc_size;
        
for (int i=0; i<chunksize[taskid]; i++)
    {
 if (i==0) acc_chunk.push_back(acc[0]);
 else  
    acc_chunk.push_back(acc[i]-acc[i-1]);
    std::cout<<"PROC["<<taskid<<"]"<<" chunk "<<i<<" "<<"accepted= "<<acc_chunk[i]<<std::endl;
    }

        //std::cout<<" Proc "<<taskid<<" is HERE 2"<<std::endl;
		int totNacc = 0;
		#ifdef PARALLEL
		MPI_Reduce(&Nacc[taskid], &totNacc, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
		#endif
        
//      for (taskid=0; taskid<numtasks; taskid++) {
//		totNacc += Nacc[taskid];
std::cout<<"MPI_Reduce PROC["<<taskid<<"] Nacc= "<<Nacc[taskid]<<" totNacc= "<<totNacc<<std::endl;
//		}
        #ifdef PARALLEL
        MPI_Bcast(&totNacc, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        #endif
        
        std::cout<<"PROC["<<taskid<<"] believes that totNacc is "<<totNacc<<std::endl;

        // Each processor passes its number of accepted nuclei to the all other processor (will be used to generate the tr1::array to pass results to the master)
        #ifdef PARALLEL
        //MPI_Allgather(&Nacc[taskid], 1, MPI_INT, &Nacc[taskid], 1, MPI_INT, MPI_COMM_WORLD);
     
        if (taskid > MASTER){
            int dest = MASTER;
            MPI_Send(&Nacc[taskid], 1, MPI_INT, dest, tag2, MPI_COMM_WORLD);
        }
        if (taskid == MASTER) {
            for (int i=1; i<numtasks; i++) {
                int source = i;
                MPI_Recv(&Nacc[i], 1, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
            }
        }
        
        MPI_Bcast(&Nacc, numtasks, MPI_INT, MASTER, MPI_COMM_WORLD);
       
        for (int j = 0; j<numtasks; j++) {
            std::cout<<"Proc ["<<taskid<<"] believes that Proc ["<<j<<"] accepted "<<Nacc[j]<<" tasks"<<std::endl;
        }
        #endif
        
printf(" Ktstep=%d, totNacc= %d, number moves accepted by P%d : %d\n", Ktstep, totNacc, taskid, acc_size);	
        // each processor creates an tr1::array as long as the total number of accepted particles (using MPI broadcast and/or reduce sum) and populates its part
        // quantitites to be passed are:
		// (1) -> x    (2) -> y     (3) -> z    (4)-> diam    (5)-> ptype    (6) -> shape    (7) -> o1     (8) -> o2     (9) -> o3
		int nqtt = 9;   //see numbers on line just above this one..
        double data_acc[totNacc*nqtt];
        int addoffset = 0;
        
        for (int i=0; i<taskid; i++)
        {
            addoffset += Nacc[i];
        }
std::cout<<"Addoffset on proc["<<taskid<<"] is "<<addoffset<<" Nacc["<<taskid<<"]= "<<Nacc[taskid]<<" totNacc= "<<totNacc<<std::endl;
       
        int intoff = 0;
		for (int i = 0; i<chunksize[taskid];i++) {
		    for (int j=0; j<acc_chunk[i]; j++){
                int pos = (addoffset + j) * nqtt;
                data_acc[pos]	=	x_acc[intoff+j];
                data_acc[pos+1]	=	y_acc[intoff+j];
                data_acc[pos+2]	=	z_acc[intoff+j];
                data_acc[pos+3]	=	diam_acc[intoff+j];
                data_acc[pos+4]	=	(double)ptype_acc[intoff+j];
                data_acc[pos+5]	=	(double)shape_acc[intoff+j];
                data_acc[pos+6]	=	o1_acc[intoff+j];
                data_acc[pos+7]	=	o2_acc[intoff+j];
                data_acc[pos+8]	=	o3_acc[intoff+j];
std::cout<<"Data_acc PROC["<<taskid<<"] chunk "<<i<<" ptype_acc["<<i<<"]= "<<data_acc[pos+4]<<" x_acc="<<data_acc[pos]<<" diam_acc= "<<data_acc[pos+3]<<" totNacc= "<<totNacc<<" acc_chunk= "<<acc_chunk[i]<<" Addoffset= "<<addoffset<<" intoff= "<<intoff<<" pos= "<<pos<<std::endl;
		    }
            addoffset += acc_chunk[i];
            intoff += acc_chunk[i];
		}
		// the master collects all the chunks of data_acc and compiles them into new particles
        int poff[numtasks];
        for (int i=0; i<numtasks; i++) {poff[i] = 0;}

        for (int i=0; i<taskid; i++) {poff[taskid] += Nacc[i];}
        
#ifdef PARALLEL
        if (taskid > MASTER){
            int dest = MASTER;
            MPI_Send(&poff[taskid], 1, MPI_INT, dest, tag2, MPI_COMM_WORLD);
           // std::cout<<"MPI_Send data_acc from PROC["<<taskid<<"] to Master offset= "<<addoffset<<" chunksize= "<<chunksize[taskid]<<std::endl;
        }
        if (taskid == MASTER) {
            for (int i=1; i<numtasks; i++) {
                int source = i;
                MPI_Recv(&poff[i], 1, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                //std::cout<<"MPI_Recv data_acc from PROC["<<i<<"] offset= "<<addoffset<<" chunksize= "<<chunksize[i]<<std::endl;
            }
        }
#endif
        
        #ifdef PARALLEL
        if (taskid > MASTER){
            int dest = MASTER;
            MPI_Send(&data_acc[poff[taskid]*nqtt], Nacc[taskid]*nqtt, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
std::cout<<"MPI_Send data_acc from PROC["<<taskid<<"] to Master offset= "<<poff[taskid]<<" chunksize= "<<chunksize[taskid]<<std::endl;
        }
        if (taskid == MASTER) {
            for (int i=1; i<numtasks; i++) {
                int source = i;
                MPI_Recv(&data_acc[poff[i]*nqtt], Nacc[i]*nqtt, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
std::cout<<"MPI_Recv data_acc from PROC["<<i<<"] offset= "<<poff[i]<<" chunksize= "<<chunksize[i]<<std::endl;
            }
        }
		#endif

        if(taskid==MASTER){
            for (int pos=0; pos<totNacc*nqtt; pos=pos+nqtt){

                std::cout<<"Data_acc after SEND-REC PROC["<<taskid<<"] chunk "<<taskid<<" ptype_acc= "<<data_acc[pos+4]<<" x_acc="<<data_acc[pos]<<" diam_acc= "<<data_acc[pos+3]<<" totNacc= "<<totNacc<<" acc_chunk= "<<acc_chunk[taskid]<<" pos= "<<pos<<std::endl;
            }
        }
        
		#ifdef PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);	 //wait for all processors to have finished with their part
		#endif
        
		if (taskid == MASTER){
            
                Particle temp_part;
//            for(int k=0; k<chunksize[taskid]; k++){
                for (int i=0; i<totNacc; i++) {
                    
                temp_part.x = data_acc[i*nqtt];
				temp_part.y = data_acc[i*nqtt+1];
				temp_part.z = data_acc[i*nqtt+2];
				temp_part.diam = data_acc[i*nqtt+3];
				temp_part.type = (int)data_acc[i*nqtt+4];
				temp_part.dens = parttype[temp_part.type].dens;
				temp_part.mass = temp_part.dens * M_PI * temp_part.diam * temp_part.diam * temp_part.diam / 6.; //yg = 10-24 g
				int sh_id = (int)data_acc[i*nqtt+5];
				temp_part.shape = parttype[temp_part.type].nuclshape[sh_id];
				temp_part.n_attr = parttype[temp_part.type].n_attr[sh_id];
                    

				for (int j=0; j<MAX_n_shape_attr; j++) {
					temp_part.sh[j]= parttype[temp_part.type].shape_attr[sh_id][j];
				}
				temp_part.o1=data_acc[i*nqtt+6];
				temp_part.o2=data_acc[i*nqtt+7];
				temp_part.o3=data_acc[i*nqtt+8];
				temp_part.vx =0.;
				temp_part.vy=0.;
				temp_part.vz=0.;
                    
std::cout<<"Data_acc_MASTER chunk "<<" ptype= "<<temp_part.type<<" x= "<<temp_part.x<<" diam= "<<temp_part.diam<<" dens= "<<temp_part.dens<<" nattr= "<<temp_part.n_attr<<std::endl;
				parts.push_back(temp_part);
                    
//            }
			}
		}
        
		// the master assigns the correct interaction type to each new particle (based on p_type, diam, internal chemistry measure)
        if (taskid == MASTER){
		}
		
        // the master writes a new input config file for LAMMPS
		// Igor, you should already have this done
        if (taskid == MASTER){
		}
		
		// if needed, the master writes also an xyz and dump file as results of the kinetic simulation
        if (taskid == MASTER){
		}
		
		// the master runs lammps to minimize potential energy or to equilibrate in NVT, or NVE, or NPT
		if (taskid == MASTER){
		}
    
        // if needed, the master writes also an xyz and dump file as results of the kinetic simulation
        if (taskid == MASTER){
		}


		// COMPLETE DISSOLUTION 
		
		// The existing particles may disappear completely at once, with a process that is qualitatively inverse to nucleation
		// At this point, we should have a LAMMPS file with all energies per Atom2. We can use the negative of them to construct
		//	a criterio similar to the one for nucleation to get a characteristic dissolution time and choose which particles to
		//	erase
		// Once we have the nucleation and this step implemented, we are good to start grinding results!
		
		// GROWTH + DISSOLUTION  (PRELIMINARY DRAFT... WE WILL NEED TO RETHINK THIS LATER ON...)
		
		// start loop of N time steps
		
		// create a copy of all type 2 particles with slightly increased diameter and add it to the
		// xyz input file for lammps giving them type 3
		// Define a lammps tabular potential where the energy of interaction is zero
		// zero if the interparticle distance is zero.
		// Enable interactions 1-4 ; 4-4 ; 4-3 ; 1-3 (this latter will exclude the particle with its non-grown self
		//      thanks to the tabular potential
		// prints the energy per Atom2 (sum of all interactions per Atom2) to a file
		
		
		// in PARALLEL (AT LOWER LEVEL):
		// read the energy-per-Atom2 and xyz files and extract only the energies for type 4 particles
		//      These are effectively Delta-interaction-energies due to a single particle small growth
		// compute time associated to the attributed growth (via number of chemical reactions per unit surface
		// requiired. Divide growth by time and get growth rate for each particle
		
		// in PARALLEL (AT HIGHER LEVEL), do the same for dissolution
		
		// in PARALLEL (AT ALL LEVELS) integrate growth + dissolution rate
		
		// run lammps to minimize the total interaction energy
		
		// endo of GROWTH + DISSOLUTION loop
		
		
		// adjust solution chemistry composition (diffusion rate? constant composition?)
		
		// every N steps record the xyz config of type 1 and 2 attoms, with diameters, and dump quantitities
    

    }    // end of MAIN loop

	#ifdef PARALLEL
	MPI_Finalize();
	#endif

*/
    
}   /* end of main */


// A FUNCTION TO READ THE attoms IN SIM, MOLECULES IN SOLUTION, AND POSSIBLE CHEMICAL REACTIONS 
void read_chem(std::vector<Atom2> &attoms,std::vector<MolSol> &molsol, std::vector<PartType> &parttype, std::vector<ChemRX> &chemRX){
    // this function should read a file containing info regarding the overall attoms that can be present, the composition of the solution, and particles that precipitate/dissolve.
    
	// For now, fast and dirty... I will specify only a cement-like solution with only one CSH and CH as possible particle..
	
	
    //  All attoms
    //  - number of Atom2 species in the whole simulation (e.g. 4 : Ca, Si, H, O in a C3S+H20 system)
    //  - name of each Atom2 species
    
	Atom2 temp_Atom2;
	temp_Atom2.name = "Ca";
	temp_Atom2.mass = 40.078 / AN * 1.E+24;  //in yg, i.e. 10^-24 grams 
	attoms.push_back(temp_Atom2);
	temp_Atom2.name.clear();
	temp_Atom2.name = "Si";
	temp_Atom2.mass = 28.08550 / AN * 1.E+24;  //in yg, i.e. 10^-24 grams 
	attoms.push_back(temp_Atom2);
	temp_Atom2.name.clear();
	temp_Atom2.name = "O";
	temp_Atom2.mass = 15.9994 / AN * 1.E+24;  //in yg, i.e. 10^-24 grams 
	attoms.push_back(temp_Atom2);
	temp_Atom2.name.clear();
	temp_Atom2.name = "H";
	temp_Atom2.mass = 1.00794 / AN * 1.E+24;  //in yg, i.e. 10^-24 grams 
	attoms.push_back(temp_Atom2);
	temp_Atom2.name.clear();
	
    //  Solution
    //  - number of molecular species that can be in solution (e.g. 4: Ca, H2SiO4 , H2O, OH)
    //  - name of each molecular species
    //  - composition of each molecular species in solution, i.e. number of attoms of each possible species composing the molecule
    //  - attributes of the molecules in solution: all what is needed to compute the free energy of the solution
    //  - initial concentration of each molecular species in solution (number per unit volume)
    
	MolSol temp_mol;
	temp_mol.name="Ca2+";
	temp_mol.atinm.push_back(1);
	temp_mol.atinm.push_back(0);
	temp_mol.atinm.push_back(0);
	temp_mol.atinm.push_back(0);
	temp_mol.conc = 20.E-3 * (1.E-24 * AN); //molecules per nm3... the value in moles times E-24 to convert per-liters to per-nm^3, and AN is the avogadro's number
	temp_mol.Ei = 0.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.Es = 100.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.q = 2.;
	molsol.push_back(temp_mol);
	temp_mol.name.clear();
	temp_mol.atinm.clear();
	
	temp_mol.name="H2SiO4";
	temp_mol.atinm.push_back(0);
	temp_mol.atinm.push_back(1);
	temp_mol.atinm.push_back(4);
	temp_mol.atinm.push_back(2);
	temp_mol.conc = 20.E-6 * (1.E-24 * AN);
	temp_mol.Ei = 600.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.Es = 100.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.q = -2.;
	molsol.push_back(temp_mol);
	temp_mol.name.clear();
	temp_mol.atinm.clear();
	
	temp_mol.name="H2O";
	temp_mol.atinm.push_back(0);
	temp_mol.atinm.push_back(0);
	temp_mol.atinm.push_back(1);
	temp_mol.atinm.push_back(2);
	temp_mol.conc = 55.56 * (1.E-24 * AN);
	temp_mol.Ei = 200.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.Es = 100.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.q = 0.;
	molsol.push_back(temp_mol);
	temp_mol.name.clear();
	temp_mol.atinm.clear();
	
	temp_mol.name="OH";
	temp_mol.atinm.push_back(0);
	temp_mol.atinm.push_back(0);
	temp_mol.atinm.push_back(1);
	temp_mol.atinm.push_back(1);
	temp_mol.conc = (40.E-3 -40.E-6) * (1.E-24 * AN);
	temp_mol.Ei = 100.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.Es = 100.;	// yg * nm^2/ps^2 =  10-27 Kg * 10 -18 m^2 / 10-24 s^2 = 10^-21 Joules = 0.00624150934 eV
	temp_mol.q = -1.;
	molsol.push_back(temp_mol);
	temp_mol.name.clear();
	temp_mol.atinm.clear();
	
    //  Particle types
    //  - number of possible particles types that can precipitate (e.g 3: CaOH2 , CSH(II), Substrate)
    //  - name of each particle type
    //  - attributes of each particle type:
    //      . number of attoms of each species per unit volume of particle
    //      . internal energy per unit volume (could be reassigned to each particle during simulation)
    //      . something about density of states and/or vibration freqs,
    //      . surface energy (in future this should depend on a solution chemistry that might change during the simulation. That could require in future a separate function to compute them..)
    //      . list of active chemical reactions for: (1) nucleation (2) growth (3) dissolution (4) leaching (5) sequestration (6) solid-solid mass transfer
    //  - number of possible nuclei for each particle type
    //  - size, shape, and orientation of possible nuclei for each particle type
	
	std::array<double, 3>	tarr;
	std::array<double, 6>	tarr2;
	
	PartType temp_ptype;
	temp_ptype.name = "CSH(II)";
	temp_ptype.ptid = 1;    //whatever you please. This is NOT the id that will be passed to LAMMPS.
	temp_ptype.atins.push_back(2);
	temp_ptype.atins.push_back(1);
	temp_ptype.atins.push_back(9);
	temp_ptype.atins.push_back(10);
	temp_ptype.volmol = (161.2 * 1.E+21 / AN);  //from Bullard, JAmCerSoc 2008, 161.2 cm3/mol
	temp_ptype.Eisp = 1.E+6; // 10^-21 J/nm3 = 0.00624150934 eV/nm^3. Invented value here
	temp_ptype.Essp = 1.E+4; // 10^-21 J/nm3 = 0.00624150934 eV/nm^2. Invented value here
	temp_ptype.dens = 2650.; // Kg m-3 = 10-24 g nm-3 = yg nm-3
	temp_ptype.RXnuc.push_back(0);	//the pushed-back id refers to the chemical reactions specified below
	temp_ptype.RXgrow.push_back(0);
	temp_ptype.RXdiss.push_back(1);
	// in this example, I will try nuclei with diameter 0.9 nm and 1.6 nm
	temp_ptype.nuclsize.push_back(0.9);
	//temp_ptype.nuclsize.push_back(1.6);
	// in this example I will allow two possible orientations of the nuclei, totally invented
	tarr[0]=0.1; tarr[1]=0.4; tarr[2]=0.87;		
	temp_ptype.nuclori.push_back( tarr );
	//tarr[0]=-0.5; tarr[1]=0.34; tarr[2]=0.7;
	//temp_ptype.nuclori.push_back( tarr );
	// in this example I will allow two possible shapes of the nuclei, totally invented
	temp_ptype.nuclshape.push_back("Sphere");
	tarr2[0]=1.; tarr2[1]=0.; tarr2[2]=0.; tarr2[3]=0.; tarr2[4]=0.; tarr2[5]=0.;	
	temp_ptype.shape_attr.push_back( tarr2 );  //this should be managed more smartly when reading an input file..
	temp_ptype.n_attr.push_back(1);
	//temp_ptype.nuclshape.push_back("Ellipsoid");
	//tarr2[0]=1.; tarr2[1]=2.; tarr2[2]=0.5; tarr2[3]=0.; tarr2[4]=0.; tarr2[5]=0.;
	//temp_ptype.shape_attr.push_back( tarr2 );
	//temp_ptype.n_attr.push_back(3);
	parttype.push_back(temp_ptype);
	temp_ptype.name.clear();
	temp_ptype.atins.clear();
	temp_ptype.RXnuc.clear();
	temp_ptype.RXgrow.clear();
	temp_ptype.RXdiss.clear();
	temp_ptype.nuclsize.clear();
	temp_ptype.nuclori.clear();
	temp_ptype.nuclshape.clear();
	temp_ptype.shape_attr.clear();
	temp_ptype.n_attr.clear();
	
    /*
	temp_ptype.name = "CH";
	temp_ptype.ptid = 8;     //whatever you please. This is NOT the id that will be passed to LAMMPS.
	temp_ptype.atins.push_back(1);
	temp_ptype.atins.push_back(0);
	temp_ptype.atins.push_back(2);
	temp_ptype.atins.push_back(2);
	temp_ptype.volmol = (33.08 * 1.E+21 / AN);  //from Bullard, JAmCerSoc 2008, 33.08 cm3/mol 
	temp_ptype.Eisp = 1.E+6; // 10^-21 J/nm3 = 0.00624150934 eV/nm^3. Invented value here
	temp_ptype.Essp = 1.E+4; // 10^-21 J/nm3 = 0.00624150934 eV/nm^2. Invented value here
	temp_ptype.dens = 2260.; // Kg m-3 = 10-24 g nm-3 = yg nm-3
	temp_ptype.RXnuc.push_back(2);	//the pushed-back id refers to the chemical reactions specified below
	temp_ptype.RXgrow.push_back(2);
	temp_ptype.RXdiss.push_back(3);
	// in this example, I will try nuclei with diameter 0.9 nm and 1.6 nm
	temp_ptype.nuclsize.push_back(1.1);
	temp_ptype.nuclsize.push_back(0.6);
	// in this example I will allow one possible orientation of the nuclei, totally invented
	tarr[0]=0.2; tarr[1]=-0.4; tarr[2]=0.87;		
	temp_ptype.nuclori.push_back( tarr );
	// in this example I will allow one possible shape of the nuclei: spherical
	temp_ptype.nuclshape.push_back("Sphere");
	tarr2[0]=1.; tarr2[1]=0.; tarr2[2]=0.; tarr2[3]=0.; tarr2[4]=0.; tarr2[5]=0.;	
	temp_ptype.shape_attr.push_back( tarr2 );  //this should be managed more smartly when reading an input file..
	temp_ptype.n_attr.push_back(1);
	parttype.push_back(temp_ptype);
	temp_ptype.name.clear();
	temp_ptype.atins.clear();
	temp_ptype.RXnuc.clear();
	temp_ptype.RXgrow.clear();
	temp_ptype.RXdiss.clear();
	temp_ptype.nuclsize.clear();
	temp_ptype.nuclori.clear();
	temp_ptype.nuclshape.clear();
	temp_ptype.shape_attr.clear();
	temp_ptype.n_attr.clear();
    */
	
    
    /*
	temp_ptype.name = "substrate";
	temp_ptype.ptid = 4;     //whatever you please. This is NOT the id that will be passed to LAMMPS.
	temp_ptype.atins.push_back(0);
	temp_ptype.atins.push_back(0);
	temp_ptype.atins.push_back(0);
	temp_ptype.atins.push_back(0);
	temp_ptype.volmol = (100 * 1.E+21 / AN);   //invented..
	temp_ptype.Eisp = 0.E+6; // 10^-21 J/nm3 = 0.00624150934 eV/nm^3. Invented value here
	temp_ptype.Essp = 0.E+4; // 10^-21 J/nm3 = 0.00624150934 eV/nm^2. Invented value here
	temp_ptype.dens = 1000.; // Kg m-3 = 10-24 g nm-3 = yg nm-3
	// no reactions allowed for the substrate, hence no RXnuc, RXgrow, RXdiss to be specified
	/*temp_ptype.RXnuc.push_back(2);	//the pushed-back id refers to the chemical reactions specified below
	temp_ptype.RXgrow.push_back(2);
	temp_ptype.RXdiss.push_back(3);*/
	// no nucleation allowed, so nothing about the nuclei must be specified
	/*temp_ptype.nuclsize.push_back(0.9);
	temp_ptype.nuclsize.push_back(1.6);
	// neither orientation of the nuclei
	tarr[0]=0.2; tarr[1]=-0.4; tarr[2]=0.87;		
	temp_ptype.nuclori.push_back( tarr );
	// in this example I will allow one possible shape of the nuclei: spherical
	temp_ptype.nuclshape.push_back("Sphere");
	tarr2[0]=1.; tarr2[1]=0.; tarr2[2]=0.; tarr2[3]=0.; tarr2[4]=0.; tarr2[5]=0.;	
	temp_ptype.shape_attr.push_back( tarr2 );  //this should be managed more smartly when reading an input file..
	temp_ptype.n_attr.push_back(1);*/
	/*
    parttype.push_back(temp_ptype);
	temp_ptype.name.clear();
	temp_ptype.atins.clear();
	temp_ptype.RXnuc.clear();
	temp_ptype.RXgrow.clear();
	temp_ptype.RXdiss.clear();
	temp_ptype.nuclsize.clear();
	temp_ptype.nuclori.clear();
	temp_ptype.nuclshape.clear();
	temp_ptype.shape_attr.clear();
	temp_ptype.n_attr.clear();
     */
	
    //  Chemical reactions
    //  invlving molecules in solutions and attoms in the solid (can be extended to include also solid-solid mass transfers)
    //  - number of possible chemical reactions
    //  - stoichiometric coefficients (this should say how many molecules are going in solution and how many attoms are produced in the solid particle)
    //  - activation energy of the reaction in isolated conditions
	//  - change of volume of solid induced by one reaction
    
	ChemRX temp_chemrx;
	temp_chemrx.name="CSH(II)prec"; //all from Bullard JAmCerSoc 2008
	temp_chemrx.mtosol.push_back(-2);
	temp_chemrx.mtosol.push_back(-1);
	temp_chemrx.mtosol.push_back(-3);
	temp_chemrx.mtosol.push_back(-2);
	temp_chemrx.atosld.push_back(2);
	temp_chemrx.atosld.push_back(1);
	temp_chemrx.atosld.push_back(9);
	temp_chemrx.atosld.push_back(10);
	temp_chemrx.DEiso = 1000. * (10.E+21 / AN);	//units of 0.00624150934 eV. The value in J/mol is converted to 10^-21 J by mutiplying times 10^21 and dividing by the avogadro's number 
	temp_chemrx.DVsld = 0.27;  //nm3... this is close to the molecular volume of CSH(II)
	chemRX.push_back(temp_chemrx);
	temp_chemrx.name.clear();
	temp_chemrx.mtosol.clear();
	temp_chemrx.atosld.clear();
	
    /*
	temp_chemrx.name="CSH(II)diss"; //all from Bullard JAmCerSoc 2008
	temp_chemrx.mtosol.push_back(2);
	temp_chemrx.mtosol.push_back(1);
	temp_chemrx.mtosol.push_back(3);
	temp_chemrx.mtosol.push_back(2);
	temp_chemrx.atosld.push_back(-2);
	temp_chemrx.atosld.push_back(-1);
	temp_chemrx.atosld.push_back(-9);
	temp_chemrx.atosld.push_back(-10);
	temp_chemrx.DEiso = 1000. * (10.E+21 / AN);	//units of 0.00624150934 eV. Invented value in J/mol is converted to 10^-21 J by mutiplying times 10^21 and dividing by the avogadro's number 
	temp_chemrx.DVsld = -0.27;  //nm3
	chemRX.push_back(temp_chemrx);
	temp_chemrx.name.clear();
	temp_chemrx.mtosol.clear();
	temp_chemrx.atosld.clear();
	
	
	temp_chemrx.name="CHprec"; //all from Bullard JAmCerSoc 2008
	temp_chemrx.mtosol.push_back(-1);
	temp_chemrx.mtosol.push_back(0);
	temp_chemrx.mtosol.push_back(0);
	temp_chemrx.mtosol.push_back(-2);
	temp_chemrx.atosld.push_back(1);
	temp_chemrx.atosld.push_back(0);
	temp_chemrx.atosld.push_back(2);
	temp_chemrx.atosld.push_back(2);
	temp_chemrx.DEiso = 1000. * (10.E+21 / AN);	//units of 0.00624150934 eV. Invented value in J/mol is converted to 10^-21 J by mutiplying times 10^21 and dividing by the avogadro's number 
	temp_chemrx.DVsld = 0.055;  //nm3... this is close to the molecular volume of CH
	chemRX.push_back(temp_chemrx);
	temp_chemrx.name.clear();
	temp_chemrx.mtosol.clear();
	temp_chemrx.atosld.clear();
	
	
	temp_chemrx.name="CHdiss"; //all from Bullard JAmCerSoc 2008
	temp_chemrx.mtosol.push_back(1);
	temp_chemrx.mtosol.push_back(0);
	temp_chemrx.mtosol.push_back(0);
	temp_chemrx.mtosol.push_back(2);
	temp_chemrx.atosld.push_back(-1);
	temp_chemrx.atosld.push_back(0);
	temp_chemrx.atosld.push_back(-2);
	temp_chemrx.atosld.push_back(-2);
	temp_chemrx.DEiso = 1000. * (10.E+21 / AN);	//units of 0.00624150934 eV. Invented value in J/mol is converted to 10^-21 J by mutiplying times 10^21 and dividing by the avogadro's number 
	temp_chemrx.DVsld = -0.055;  //nm3
	chemRX.push_back(temp_chemrx);
	temp_chemrx.name.clear();
	temp_chemrx.mtosol.clear();
	temp_chemrx.atosld.clear();
     */
}


// A FUNCTION TO READ THE INITIAL LAMMPS-LIKE CONFIGURATION
void read_conf(std::string c0fname , double box[], std::vector<Particle> &parts){
    // this function should read a lammps configuration including:
    //  - number of particles
    //  - number of particle types
    //  - box sizes
    //  - type of pairs to define subsequent coefficients
    //  - list of pair coefficients for each Atom2 type (WHAT IS THIS??).
    //  - type of attoms style (sphere..)
    //  - rows with label, type, diameter, density, x , y ,z , bs, bs ,bs
    
    // For now, I just write the conf here myself: fast and dirty..
    
    // box sizes (here in nanometers)
    box[0]=100;
    box[1]=100;
    box[2]=100;
    
    // number of attoms in the substrate layer
    int npart;
    double DAtom2 = 5.; //nm
    int npart_x = (int)(box[0]/DAtom2);
    int npart_y = (int)(box[1]/DAtom2);
    npart = npart_x * npart_y;
	
	//interparticle distances in the layer
	double dx = box[0]/(double)npart_x;
	double dy = box[1]/(double)npart_y;
    

	Particle temp_part;
	//create initial substrate particles
    for (int i=0 ; i<npart ; i++){
		temp_part.type = 4;	//"4" is the id that I assigned to the substrate particles beforehead
        temp_part.diam = 5.; //nm
		temp_part.dens = 2650.; // Kg m-3 = 10-24 g nm-3 = yg nm-3
		temp_part.mass = temp_part.dens * M_PI * temp_part.diam * temp_part.diam * temp_part.diam / 6.; //yg = 10-24 g
		temp_part.x = ((double)(i % npart_x) + 0.5) * dx;
		temp_part.y = ((double)( (int)(i/npart_x)) + 0.5) * dy;
		temp_part.z = 0.;
		temp_part.vx = 0.;
		temp_part.vy=0.;
		temp_part.vz=0;
		temp_part.o1=0.;
		temp_part.o2=0.;
		temp_part.o3=0.;
		temp_part.shape = "Sphere";
		temp_part.n_attr = 1;
		temp_part.sh[0]= 1.;
		temp_part.sh[1]= 0.; temp_part.sh[2]= 0.; temp_part.sh[3]= 0.; temp_part.sh[4]= 0.; temp_part.sh[5]= 0.;
		parts.push_back(temp_part);
    }
	//create 10 initial CSH(II) particles
	for (int i=0 ; i<10 ; i++){
		temp_part.type = 1;	//"1" is the id that I assigned to the CSH(II) particles beforehead
        temp_part.diam = 5.; //nm
		temp_part.dens = 2650.; // Kg m-3 = 10-24 g nm-3 = yg nm-3
		temp_part.mass = temp_part.dens * M_PI * temp_part.diam * temp_part.diam * temp_part.diam / 6.; //yg = 10-24 g
		temp_part.x = 0.01+box[0]/10.*sqrt((double)i);
		temp_part.y = 0.02+box[0]/15.*(double)i;
		temp_part.z = 20. + (double)(i%2)*30.;
		temp_part.vx = 0.;
		temp_part.vy=0.;
		temp_part.vz=0;
		temp_part.o1=0.;
		temp_part.o2=0.;
		temp_part.o3=0.;
		temp_part.shape = "Sphere";
		temp_part.n_attr = 1;
		temp_part.sh[0]= 1.;
		temp_part.sh[1]= 0.; temp_part.sh[2]= 0.; temp_part.sh[3]= 0.; temp_part.sh[4]= 0.; temp_part.sh[5]= 0.;
		parts.push_back(temp_part);
    }
	//create 2 initial CH particles
	for (int i=0 ; i<2 ; i++){
		temp_part.type = 8;	//"8" is the id that I assigned to the CH particles beforehead
        temp_part.diam = 7.; //nm
		temp_part.dens = 2260.; // Kg m-3 = 10-24 g nm-3 = yg nm-3
		temp_part.mass = temp_part.dens * M_PI * temp_part.diam * temp_part.diam * temp_part.diam / 6.; //yg = 10-24 g
		temp_part.x = 10.+ (double)i*20;
		temp_part.y = 10.+ (double)i*20;
		temp_part.z = 30. + (double)(i%2)*30.;
		temp_part.vx = 0.;
		temp_part.vy=0.;
		temp_part.vz=0;
		temp_part.o1=0.;
		temp_part.o2=0.;
		temp_part.o3=0.;
		temp_part.shape = "Sphere";
		temp_part.n_attr = 1;
		temp_part.sh[0]= 1.;
		temp_part.sh[1]= 0.; temp_part.sh[2]= 0.; temp_part.sh[3]= 0.; temp_part.sh[4]= 0.; temp_part.sh[5]= 0.;
		parts.push_back(temp_part);
    }
}


// A FUNCTION TO READ THE CHEMICAL INFO OF THE PARTICLES IN THE INITIAL CONFIGURATION (extra info to what lammps handles)
void read_Pchem(std::vector<Particle> &parts, std::vector<PartType> &parttype, std::vector<Atom2> &attoms){

	// for now, I just read the particle chemical info referring to the corresponding particle type defined previously
	// in future, for simulation restart purposes, it will be better to read also the chemical info of each particles from a file
	// saved along with the xyz coordinates. In fact, in future, the specific chemical composition of single particles might depart from
	// that in PartType as a consequence of leaching, enrichment, possibly growth and dissolution too.

	
	
	for (int i=0 ; i<parts.size() ; i++){
		
		//VOLUME AND SURFACE ARE COMPUTED ASSUMING SPHERICAL PARTICLES!!
		parts[i].vol = M_PI / 6. * parts[i].diam * parts[i].diam * parts[i].diam;
		parts[i].surf = M_PI * parts[i].diam * parts[i].diam;
		
		//number of attoms of different species
		for (int j=0; j<attoms.size(); j++ ) {
			parts[i].atinp.push_back( parts[i].vol/parttype[parts[i].typepos].volmol  * parttype[parts[i].typepos].atins[j] );
		}
		
		parts[i].Ei = parttype[parts[i].typepos].Eisp * parts[i].vol;
		parts[i].Esurf = parttype[parts[i].typepos].Essp * parts[i].surf;
	}
}



// A FUNCTION TO READ THE LIST OF POSSIBLE INTERACTIONS AND TO ASSIGN THEM TO PARTICLES TYPE AND TO CREATE SUBGROUPS OF PARTICLES FOR LAMMPS
/*void create_int_table(){
    
    // first, this should read through a file listing particle types (id or name) and corresponding interaction type
    // for example, for the case here, something like:
    //  CSH(II) - substrate : Lennard-Jones eps_value   sigma_value     (see LAMMPS)
    //  CH - CH        : Mei           eps_value   alpha_value     exp1_value  exp2_value       (see LAMMPS)
    //  CSH(II) - CSH(II)   : PolyMei       eps/sigma3_value            exp1_value  exp2_value    (currently not exising in LAMMPS)
    //  CH - substrate      : PairTab        (see LAMMPS)
    //                          dist1   energy1     force1
    //                          dist2   energy2     force2
    //                          .......
    //  CSH(II) - CH        : MDCTab      Pair Tab
    //                          dist1   energy1     force1
    //                          dist2   energy2     force2
    //                          .......
    //                        pt1_diam1    pt2_diam1     pt1_ChemInd1    pt2_ChemInd1   pt1_Mod1  pt2_Mod1  Escale1   pt1_ChemInd2    pt2ChemInd2  pt1_Mod2 pt2_Mod2
    //                        pt1_ChemInd = ratio_attoms  Ca   Si
    //                        pt2_ChemInd = rtaio_attoms  Ca   O         (more realistic, in future, could be carbonation depth or so)
    //                        pt1_Dmin  pt1_Dmax  10  pt1_ChemIndmin    pt1_ChemIndmax 5
    //                        pt2_Dmin  pt2_Dmax  4  pt2_ChemIndmin    pt2_ChemIndmax  1
    //
    //
    //  NB: MDC stands for multi-diameter-chemistry-table. It requires as inputs:
    //      - a potential type that can be LJ or Mei or PairTab with associated attributes (see above)
    //      - diameters a pair of particle_types 1 and 2 for which the energy scaling parameter Escale is known at two given
    //          values of chemical composition (chemical indicator). Escale = 1 corresponds to potential type input above.
    //      - this inputs allow the code to provide scaling
    //          E = E1 / (beta1 * Davrg1^3) * (beta*Davrg^3) /  Yeff1 * Yeff
    //              where E = Escale, beta see Masoero et al SoftMatter 2014, Davrg = average diameter, Davrg1 = D avrg of pt1_diam1
    //                  and pt2_diam1, Yeff is the arithmetic average of the elastic moduli pt1_Mod and pt2_Mod. I assume that
    //                  the moduli vary linealy between two known chemical compositions ChemInd1 and ChemInd2, both for pt1 and pt2
    //      - the type of chemical indicator is specified for both particle types involved
    //      - the range of diameters for different potential tables for pt1 is specified, along with the number of bins
    //      - same for chemical composition, and for particle type pt2
    //      - This code will create a 10*5*4*1 possible MDCTab interaction types with the original PairTab rescaled according to the current diameter+composition bin. The CSH particles will have 10*5 possible natures according to the bin they pertain to. The CH particles will have 4*1 possible bins, according to the diameter-composition bin they pertain too. Particles with size smaller than the smallest size bin will be placed in the smallest size bin. Same for other bin types.
 
    //
    // Here I will set such input by hand (fast and dirty....)
    
    Interact temp_inter;
    temp_inter.ptype1 = 1;  // the id attributed to CSH(II). It will need to be converted to a ptypepos to point to the right entry in the parttype vector
    temp_inter.ptype2 = 4;  // the id of the substrate
    temp_inter.name = "Lennard-Jones"
}
*/


// A FUNCTION TO PRINT ALL THE INPUTS AND CHECK THAT ALL IS WORKING FINE
void print_input(std::vector<Atom2> &attoms,std::vector<MolSol> &molsol, std::vector<PartType> &parttype, std::vector<ChemRX> &chemRX, double box[], std::vector<Particle> &parts){
	
	
	printf("\n %i Atom2 species specified : ", (int)attoms.size());
	for (int i=0; i<attoms.size(); i++) {printf(" %s ", attoms[i].name.c_str());}
	printf("\n");
	printf("\n Their masses are : ");
	for (int i=0; i<attoms.size(); i++) {printf(" %f ", attoms[i].mass);}
	printf("\n");
	printf("\n %i species of molecules in solution : ", (int)molsol.size());
	for (int i=0; i<molsol.size(); i++) {printf(" %s ", molsol[i].name.c_str());}
	printf("\n");
	printf("\n Their chemical formulas are : ");
	for (int i=0; i<molsol.size(); i++) {
		for (int j=0;j<attoms.size(); j++) {
			if (molsol[i].atinm[j]>0) {
				printf("%1.1f%s",molsol[i].atinm[j],attoms[j].name.c_str());
			}
		}
		printf("  ");
	}
	printf("\n");
	printf("\n Their internal and solvation energies and charges are :");
	for (int i=0; i<molsol.size(); i++) {printf("\n %f %f %f", molsol[i].Ei,molsol[i].Es,molsol[i].q);}
	printf("\n");
	printf("\n %i possible chemical reactions : ", (int)chemRX.size());
	for (int i=0; i<chemRX.size(); i++) {printf(" %s ", chemRX[i].name.c_str());}
	printf("\n");
	printf("\n The reaction formulas are : \n");
	for (int i=0; i<chemRX.size(); i++) {
		for (int j=0;j<chemRX[i].mtosol.size(); j++) {
			if (chemRX[i].mtosol[j]!=0) {
				printf("%1.1f%s +",-chemRX[i].mtosol[j],molsol[j].name.c_str());
			}
		}
		printf(" --> ");
		for (int j=0;j<chemRX[i].atosld.size(); j++) {
			if (chemRX[i].atosld[j]!=0) {
				printf("%1.1f%s",chemRX[i].atosld[j],attoms[j].name.c_str());
			}
		}
		printf("\n");
	}
	printf("\n");
	printf("\n Their activation energies are :");
	for (int i=0; i<chemRX.size(); i++) {printf(" %f ", chemRX[i].DEiso);}
	printf("\n");
	printf("\n Their associated change in solid volume are :");
	for (int i=0; i<chemRX.size(); i++) {printf(" %f ", chemRX[i].DVsld);}
	printf("\n");
	printf("\n %i possible particle types : ", (int)parttype.size());
	for (int i=0; i<parttype.size(); i++) {printf(" %s ", parttype[i].name.c_str());}
	printf("\n");
	printf("\n Their chemical formulas are : ");
	for (int i=0; i<parttype.size(); i++) {
		bool any_form = false;
		for (int j=0;j<attoms.size(); j++) {
			if (parttype[i].atins[j]>0) {
				any_form=true;
				printf("%1.1f%s",parttype[i].atins[j],attoms[j].name.c_str());
			}
		}
		if (!any_form) {
			printf("none");
		}
		printf("  ");
	}
	printf("\n");
	printf("\n Their associated molecular volume and specific internal and surface energies :");
	for (int i=0; i<parttype.size(); i++) {printf("\n %f %f %f ", parttype[i].volmol, parttype[i].Eisp, parttype[i].Essp);}
	printf("\n");
	printf("\n Their possible nucleation reactions are :");
	for (int i=0; i<parttype.size(); i++) {
		printf("\n");
		for (int j=0; j<parttype[i].RXnuc.size(); j++) {
			printf("%s ", chemRX[ parttype[i].RXnuc[j] ].name.c_str());
		}
		if (parttype[i].RXnuc.size()==0) {
			printf("none");
		}
	}
	printf("\n");
	printf("\n Their possible growth reactions are :");
	for (int i=0; i<parttype.size(); i++) {
		printf("\n");
		for (int j=0; j<parttype[i].RXgrow.size(); j++) {
			printf("%s ", chemRX[ parttype[i].RXgrow[j] ].name.c_str());
		}
		if (parttype[i].RXgrow.size()==0) {
			printf("none");
		}
	}
	printf("\n");
	printf("\n Their possible dissolution reactions are :");
	for (int i=0; i<parttype.size(); i++) {
		printf("\n");
		for (int j=0; j<parttype[i].RXdiss.size(); j++) {
			printf("%s ", chemRX[ parttype[i].RXdiss[j] ].name.c_str());
		}
		if (parttype[i].RXdiss.size()==0) {
			printf("none");
		}
	}
	printf("\n");
	printf("\n Their possible nuclei sizes are :");
	for (int i=0; i<parttype.size(); i++) {
		printf("\n");
		for (int j=0; j<parttype[i].nuclsize.size(); j++) {
			printf("%f ", parttype[i].nuclsize[j]);
		}
		if (parttype[i].nuclsize.size()==0) {
			printf("none");
		}
	}
	printf("\n");
	printf("\n Their possible nuclei orientations are :");
	for (int i=0; i<parttype.size(); i++) {
		printf("\n");
		for (int j=0; j<parttype[i].nuclori.size(); j++) {
			printf("{ %f %f %f }  ", parttype[i].nuclori[j][0], parttype[i].nuclori[j][1], parttype[i].nuclori[j][2] );
		}
		if (parttype[i].nuclori.size()==0) {
			printf("none");
		}
	}
	printf("\n");
	printf("\n Their possible nuclei shapes are :");
	for (int i=0; i<parttype.size(); i++) {
		printf("\n");
		for (int j=0; j<parttype[i].nuclshape.size(); j++) {
			printf("%s : { ",parttype[i].nuclshape[j].c_str());
			for (int k=0; k<parttype[i].n_attr[j]; k++) {
				printf(" %f  ", parttype[i].shape_attr[j][k]);
			}
			printf("}    ");
		}
		if (parttype[i].nuclshape.size()==0) {
			printf("none");
		}
	}
	printf("\n");
	
	
	// all particles
	printf("\nBox sizes are %f %f %f \n",box[0],box[1],box[2]);
	printf("\n number of particles is %i \n",(int)parts.size());
	printf("\n ID  type   type_id      type_pos   ");
	for (int i=0; i<attoms.size(); i++) { printf("%s    ",attoms[i].name.c_str()); }
	printf("vol   surf    diam    mass    dens     x    y    z    vx     vy     vz     Ei    Esurf    o1     o2     o3    shape     n_attr   ");
	for (int i=0; i<MAX_n_shape_attr; i++) { printf("sh%i    ",i); }
	for (int i=0 ; i<parts.size() ; i++){
		printf("\n %i   %s   %i    %i    ",i,parttype[parts[i].typepos].name.c_str(),parts[i].type,parts[i].typepos);
		for (int j=0; j<attoms.size(); j++) { printf("%f    ",parts[i].atinp[j]); }
		printf(" %f   %f   %f    %f     %f   %f    %f     %f   ",parts[i].vol,parts[i].surf,parts[i].diam,parts[i].mass,parts[i].dens,parts[i].x,parts[i].y,parts[i].z);
		printf(" %f   %f   %f    %f     %f   %f    %f     %f   ",parts[i].vx,parts[i].vy,parts[i].vz,parts[i].Ei,parts[i].Esurf,parts[i].o1,parts[i].o2,parts[i].o3);
		printf(" %s   %i   ",parts[i].shape.c_str(),parts[i].n_attr);
		for (int j=0; j<MAX_n_shape_attr; j++) { printf("%f    ",parts[i].sh[j]); }
	}
	printf("\n");
}
