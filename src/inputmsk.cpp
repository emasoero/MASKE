#include "inputmsk.h"
#include "error.h"
#include "universe.h"
#include "DTnucleate.h"
#include "lammpsIO.h"
#include "simbox.h"
#include "particles.h"
#include "chemistry.h"
#include "interactions.h"
#include "solution.h"
#include "fix.h"
#include "krun.h"
#include "output.h"
#include "relax.h"
#include "block.h"
#include "store.h"

#include <unistd.h>   //just for the sleep() function
//#include "randm.h"

//needed to convert strings to lower case
#include <algorithm>


//for cluster Rocket...
#include <string.h>


using namespace MASKE_NS;

Inputmsk::Inputmsk(MASKE *maske, int argc, char **argv) : Pointers(maske)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating inputmsk object\n");
    
    if (argc <2) {
        std::string msg = "ERROR: no inputmsk file specified\n";
        error->errsimple(msg);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    isrestart = false;
    fname = argv[1];
    time_strict=false;
    time_loose=false;
    time_first=false;
    time_last=false;
    ru = 1.;
    foundNsteps=false;
    foundSubcomm=false;
    foundstep=false;
    foundtempo=false;

}



// ---------------------------------------------------------------
// Class destructor
Inputmsk::~Inputmsk()
{
    
}







// ---------------------------------------------------------------
// Reading inputmsk file (not from restart)
void Inputmsk::file()
{
	
    if (me == MASTER) fprintf(screen, "\nLooking for \"MSKrestart.dat\" file in current folder...\n");
    
    std::string read_string;
    
    std::ifstream inFile("MSKrestart.dat");
    
    
    
    // Is this a restart or not?
    if (!inFile.is_open()) {
        inFile.open (fname.c_str(), std::ifstream::in);
        if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" not found. Proceeding to input file \"%s\"...\n",fname.c_str());
    }
    else{
        if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" file found. Checking if it is active...\n");
        //look for the STATUS line, which must be the first non-comment line
        bool isrestartfound = false;
        while (!inFile.eof() && !isrestartfound) {
            inFile >> read_string;
            std::transform(read_string.begin(), read_string.end(), read_string.begin(), ::tolower);  //convert string lower case
            if (strncmp(read_string.c_str(), "#", 1) == 0) {
                std::getline (inFile, totrash); // lines starting with # are comments, which are disregarded here
            }
            else if (strcmp(read_string.c_str(), "status") == 0) {
                inFile >> read_string;
                std::transform(read_string.begin(), read_string.end(), read_string.begin(), ::tolower);  //convert string lower case
                if (strcmp(read_string.c_str(), "active") == 0){
                    isrestart = true;
                    isrestartfound = true;
                    if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" file is active; let's go with it...\n");
                }
                else{
                    if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" is not active (to make it active, the first non-comment line in \"MSKrestart.dat\" should be: \n STATUS ACTIVE \nProceeding to inputcprs file  \"%s\"...\n",fname.c_str());
                    inFile.close();
                    inFile.open (fname.c_str(), std::ifstream::in);
                    isrestartfound = true;
                }
            }
            else {
                if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" is not active (to make it active, the first non-comment line in \"MSKrestart.dat\" should be: \n STATUS ACTIVE \nProceeding to inputcprs file...\n");
                inFile.close();
                inFile.open (fname.c_str(), std::ifstream::in);
                isrestartfound = true;
            }
        }
    }
    
	
	//crash if the input file cannot be read
	if (!inFile.is_open()) {
        err_msg = "ERROR: cannot read file \""+fname+"\"";
        error->errsimple(err_msg);
	}
    MPI_Barrier(MPI_COMM_WORLD);
	
    
    
    
    
   // READ FILE (all processors need to know this)
    
    while (!inFile.eof()) {
        //inFile >> read_string;
        
        /*if (foundSubcomm){
            fprintf(screen,"\n I am processor %d , part of subcomm %s , and I am waiting here. Last command was %s \n",me,(universe->SCnames[universe->color]).c_str(),read_string.c_str());
        }
         */
        MPI_Barrier(MPI_COMM_WORLD);
        /*if (foundSubcomm){
            fprintf(screen,"\n I am processor %d , part of subcomm %s , I passed the barrier after last command was %s \n",me,(universe->SCnames[universe->color]).c_str(),read_string.c_str());
        }
         */
        std::getline (inFile, read_string);
        
        if (!read_string.empty()){
            std::istringstream ss(read_string);
            ss >> firstWord;
            
            // if this is not a loop or a block command, execute the line
            // if it is a loop, create a vector containing all line commands in the loop and executes them for the desired number of times
            // if it is a block, record the block lines in object and move on
            if (strcmp(firstWord.c_str(),"loop")==0) {
                // read the loop name and extremes. Record following lines in vector until the end_of_loop line is found
                std::string iname,imaxname,lmpyn,loopname;
                int i0,imax;
                ss >> loopname>>iname >>imaxname>>imax>>lmpyn;
                if (strcmp(lmpyn.c_str(),"lmp_yes")==0 && lammpsIO->lammps_active) {
                    std::ostringstream todo;
                    // todo << "variable "<<iname<<" equal "<<i0;
                    //lammpsIO->lammpsdo(todo.str());
                    //todo.str("");
                    todo << "variable "<<imaxname<<" equal "<<imax;
                    lammpsIO->lammpsdo(todo.str());
                }
                bool foundendloop = false;
                // record the operations in the loop
                std::vector<std::string> opers;
                
                while (!foundendloop) {
                    std::getline (inFile, read_string);
                    std::istringstream css(read_string);
                    std::string endloopname;
                    css >> word>>endloopname;
                    if (strcmp(word.c_str(),"end_of_loop")!=0 || strcmp(endloopname.c_str(),loopname.c_str())!=0)  opers.push_back(read_string);
                    else foundendloop=true;
                }
                // execute the loop for the wanted number of times
                for (int i=0; i<imax; i++) {
                    //if (me == MASTER) fprintf(screen," \n\n\n    EDDIG JO   \n\n\n");
                    //if (me == MASTER) fprintf(screen," \n\n\n    opers size = %d   \n\n\n",opres.soze);
                    if (strcmp(lmpyn.c_str(),"lmp_yes")==0 && lammpsIO->lammps_active) {
                        std::ostringstream todo;
                        todo << "variable "<<iname<<" equal "<<i+1;
                        lammpsIO->lammpsdo(todo.str());
                    }
                    for (int j=0; j<opers.size(); j++) {
                        //if (me == MASTER) fprintf(screen," \n\n\n    EDDIG JO %d   \n\n\n",j);
                        execline(opers[j]);
                    }
                }
            }
            else if (strcmp(firstWord.c_str(), "block") == 0){
                std::string read_string3, read_string4;
                ss>>read_string2>>read_string3;
                if (me==MASTER) {
                    fprintf(screen,"\n inputmsk -- Adding block:  %s %s\n",read_string2.c_str(),read_string3.c_str());
                }
                
                bool endblock = false;
                int nread = 0;
                if (strcmp(read_string2.c_str(), "lattice") == 0){
                    block->add_lattice(read_string3);
                    while (endblock == false && nread < 10000){
                        std::getline(inFile, read_string4);
                        if (strcmp(read_string4.c_str(), "end_block") == 0){
                            endblock = true;
                            block->LatLammps.push_back(block->tempSvec);
                        }
                        else {
                            block->add_Latline(read_string4);
                        }
                        nread++;
                    }
                }
                else if (strcmp(read_string2.c_str(), "region") == 0){
                    block->add_region(read_string3);
                    while (endblock == false && nread < 10000){
                        std::getline(inFile, read_string4);
                        if (strcmp(read_string4.c_str(), "end_block") == 0){
                            endblock = true;
                            block->RegLammps.push_back(block->tempSvec);
                        }
                        else {
                            block->add_Regline(read_string4);
                        }
                        nread++;
                    }
                }
                
                if (!endblock){
                    std::string msg = "ERROR: lattice or region block named \""+read_string3+"\" does not terminate with end_block keyword (or is >10000 lines long: to enable >10000 lines blocks change inputmsk.cpp";
                    error->errsimple(msg);
                }
                
                
                // verify how blocks are read via screen output
                if (me==MASTER){
                    fprintf(screen,"\n inputmsk -- Proc %d recorded block %s called %s with following lammps commands:\n",me,read_string2.c_str(),read_string3.c_str());
                    if (strcmp(read_string2.c_str(), "lattice")==0){
                        int last  = block->LatLammps.size();
                        int lastlast = block->LatLammps[last-1].size();
                        for (int i=0; i<lastlast; i++){
                            fprintf(screen,"%s\n",block->LatLammps[last-1][i].c_str());
                        }
                        fprintf(screen,"%s\n",block->LatMinMod[last-1].c_str());
                        fprintf(screen,"%s\n",block->LatMinCmd[last-1].c_str());
                    }
                     if (strcmp(read_string2.c_str(), "region")==0){
                         int last  = block->RegLammps.size();
                         int lastlast = block->RegLammps[last-1].size();
                         for (int i=0; i<lastlast; i++){
                             fprintf(screen,"%s\n",block->RegLammps[last-1][i].c_str());
                         }
                    }
                }
                
                if (strcmp(read_string2.c_str(), "lattice")==0){
                    if (block->MinCmdFound == false || block->MinModFound == false || block->LatFound==false) {
                        std::string msg = "ERROR: lattice block named \""+read_string3+"\" must contain a minimize, a min_modify, and a lattice commmand for LAMMPS";
                        error->errsimple(msg);
                    }
                }
                if (strcmp(read_string2.c_str(), "region")==0){
                    if (block->RegFound == false) {
                        std::string msg = "ERROR: region block named \""+read_string3+"\" must contain a region commmand for LAMMPS";
                        error->errsimple(msg);
                    }
                }
                if (block->BlockFound){
                    std::string msg = "ERROR: region block named \""+read_string3+"\" either contains another nested block (forbidden) or is not ended correctly with an end_block command";
                    error->errsimple(msg);
                }
                
            }
	    // Check for store multi command because we need to parse multiple lines
            else if (strcmp(firstWord.c_str(), "store") == 0) {
	      std::string type;
	      ss >> type;
	      if (type == "multi") {
		std::string name;
		ss >> name;
		store->MulNames.push_back(name);
		std::string triple_quote;
		ss >> triple_quote;
		if (triple_quote != "\"\"\"")
		  error->errsimple("ERROR: couldn't find triple quote in multi-line store command");
		int nline = 0;
		std::vector<std::string> cmds;
		while (true) {
		  std::string line;
		  std::getline(inFile, line);
		  if (line.empty())
		    continue;
		  if (line.rfind("\"\"\"", 0) == 0)
		    break;
		  cmds.push_back(line);
		  ++nline;
		  if (nline > 1000000)
		    error->errsimple("ERROR: too many lines in multi-line store command");
		}
		store->MulCmd.push_back(cmds);
	      }
	      else
		execline(read_string);
	    }
            else {
                execline(read_string);
            }
        }
    }
    inFile.close();
    
    
    // Check that all what needed to be foound in the inputcprs file has been found indeed
    // TO BE POLISHED
    /*if (!foundNsteps) {
        std::string msg = "ERROR: Nsteps not found in inputcprs file";
        error->errsimple(msg);
    }*/
    if (!foundSubcomm) {
        std::string msg = "ERROR: No subcommunicator found";
        error->errsimple(msg);
    }
    /*if (!foundstep && isrestart) {
        std::string msg = "ERROR: step not found in restart file";
        error->errsimple(msg);
    }
     */
    if (!foundtempo && isrestart) {
        std::string msg = "ERROR: No tempo found";
        error->errsimple(msg);
    }

    // Check that the number of processors in the subcommunicators equals the total number of processors passed to the program
    int nprocs,nps=0;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    for (int i=0; i<universe->nsc; i++) nps += universe->SCnp[i];
    if (nps != nprocs) {
        std::string msg = "ERROR: Number of processors invoked by mpiexec does not match user-defined number of processors needed to fill subcommunicators";
        error->errsimple(msg);
    }
                                      
    
    
    // Finish reading the inputcprs file
    if (me == MASTER) fprintf(screen, "\nDone reading the input (or restart) file...\n");
    MPI_Barrier(MPI_COMM_WORLD);
}




// ---------------------------------------------------------------
// reading initial configuration from xyz lammps-dump-like file
void Inputmsk::execline(std::string read_string)
{
    
    bool getout = false;
    
    // each line is read word by word. Some words are keywords after which a number of elements is expected. As soon as the symbol # is found the rest of line is trashed being a comment
    std::istringstream lss(read_string);
    while (lss >> word && !getout) {
        if (strncmp(word.c_str(), "#", 1) == 0) break;
        else if (strcmp(word.c_str(), "Nsteps") == 0) {
            lss >> msk->Nsteps;
            if (me == MASTER) fprintf(screen, "\nNsteps = %d", msk->Nsteps);
            foundNsteps = true;
        }
        else if (strcmp(word.c_str(), "step") == 0 && isrestart) {
            lss >> msk->step;
            msk->doublestep = (double) msk->step;
            if (me == MASTER) fprintf(screen, "\nlast step in restart was %d", msk->step);
            foundstep = true;
        }
        else if (strcmp(word.c_str(), "tempo") == 0 && isrestart) {
            lss >> msk->tempo;
            if (me == MASTER) fprintf(screen, "\ntempo = %f", msk->tempo);
            foundtempo = true;
        }
        else if (strcmp(word.c_str(), "subcomm") == 0) {
            lss >> universe->nsc;
            for (int i=0; i<universe->nsc; i++) {
                int np, seme;
                bool lmpyn;
                lss >> read_string >> np >> read_string2 >> seme;
                if (strcmp(read_string2.c_str(), "lmp_yes") == 0) lmpyn = 1;
                else if (strcmp(read_string2.c_str(), "lmp_no") == 0) lmpyn = 0;
                else {
                    err_msg = "ERROR: In input file, when defininf subcomm, I expected \"lmp_yes\" or \"lmp_no\", and instead I found \""+read_string2+"\"";
                    error->errsimple(err_msg);
                }
                universe->addsubcomm(read_string,np,lmpyn,seme);
                if (me == MASTER){
                    std::string temp_string;
                    temp_string = universe->SCnames[i];
                    fprintf(screen, "\nSubcommunicator \"%s\" added (with %d processors, lammps = %d, and seme = %d)",temp_string.c_str(),universe->SCnp[i],universe->SClmp[i],universe->SCseme[i]);
                }
            }
            universe->create();
            lammpsIO->create();
            foundSubcomm = true;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (strcmp(word.c_str(), "real_types") == 0) {
            int tt;
            while (lss >> tt){  msk->Rtypes.push_back(tt); }
            if (me == MASTER) {
                fprintf(screen, "Recorded real particle types: ");
                for (int i=0; i< msk->Rtypes.size(); i++) fprintf(screen, " %d ",msk->Rtypes[i]);
                fprintf(screen, "\n");
            }
        }
        else if (strcmp(word.c_str(), "trial_types") == 0) {
            int tt;
            while (lss >> tt){  msk->Ttypes.push_back(tt); }
            if (me == MASTER) {
                fprintf(screen, "Recorded trial particle types: ");
                for (int i=0; i< msk->Ttypes.size(); i++) fprintf(screen, " %d ",msk->Ttypes[i]);
                fprintf(screen, "\n");
            }
        }
        else if (strcmp(word.c_str(), "parttype") == 0){
            lss >> read_string2;
            particles->addtype(read_string2);
        }
        else if (strcmp(word.c_str(), "chemDB") == 0){
            lss >> read_string;
            chem->readDB(read_string);
        }
        else if (strcmp(word.c_str(), "kB") == 0) {
            lss >> std::scientific >> msk->kB;
            if (me == MASTER) fprintf(screen, "\n kB = %e", msk->kB);
        }
        else if (strcmp(word.c_str(), "hpl") == 0) {
            lss >> msk->hpl;
            if (me == MASTER) fprintf(screen, "\n hpl = %e", msk->hpl);
        }
        else if (strcmp(word.c_str(), "lammps") ==0)  {
            if (lammpsIO->lammps_active) {
                //fprintf(screen, "In %s LAMMMPS doing: %s\n",(universe->SCnames[universe->color]).c_str(),(lss.str()).c_str());
                //std::getline (inFile, read_string);
                //std::istringstream iss(read_string);
                std::string subcname;
                lss >> subcname;
                if (strcmp(subcname.c_str(),(universe->SCnames[universe->color]).c_str())==0 || strcasecmp(subcname.c_str(),"all")==0 ) {
                    std::string subs, newstring;
                    lss >> firstWord;
                    // add subcomm name to dump name ----
                    if (strcmp(firstWord.c_str(),"dump")==0) {
                        newstring = "dump ";
                        for (int i=0; i<4; i++) {
                            lss >> subs;
                            newstring = newstring+subs+" ";
                        }
                        lss >> subs;
                        subs = subs+"."+universe->SCnames[universe->color];
                        newstring = newstring+subs+" ";
                        while (lss >> subs) {
                            newstring = newstring+subs+" ";
                        }
                        //fprintf(screen, "\nString is \"%s\"\n",newstring.c_str());
                    }
                    else{
                        newstring = firstWord+" ";
                        while (lss >> subs) {
                            newstring = newstring+subs+" ";
                        }
                        if (strcmp(firstWord.c_str(),"thermo_style")==0) {
                            lammpsIO->lmpThSt = newstring;
                        }
                    }
                    //---
                    //fprintf(screen, "In %s LAMMMPS doing: %s\n",(universe->SCnames[universe->color]).c_str(),newstring.c_str());

                    lammpsIO->lammpsdo(newstring);
                }
                else {
                    break;   //breaks this execline function, going back to main function and reading next line
                }
            }
            else {
                getout=true; // seems to be doing the same as "break" above... not sure though.. I'll keep it as is for now
            }
        }
        else if (strcmp(word.c_str(), "simbox") == 0){
            lss>>read_string2;
            if (strcmp(read_string2.c_str(), "prism") == 0) {
                lss >> simbox->xyzlo[0] >> simbox->xyzhi[0]  >> simbox->xyzlo[1] >> simbox->xyzhi[1]  >> simbox->xyzlo[2] >> simbox->xyzhi[2]   >> simbox->xyztri[0] >> simbox->xyztri[1] >> simbox->xyztri[2];
                simbox->computeTM();
            }
            else {
                std::string msg = "ERROR: thus far only \"prism\" is accepted as a simulation box style. \""+read_string+"\" instead was found \n";
                error->errsimple(msg);
            }
        }
        else if (strcmp(word.c_str(), "pair_coeff") == 0){
            std::string read_string3;
            lss>>read_string>>read_string2;
            std::getline(lss, read_string3);
            interact->addpcoeff(read_string,read_string2,read_string3);
        }
        else if (strcmp(word.c_str(), "sol_start") == 0){
            if (me == MASTER) fprintf(screen, "\nFound sol_start\n");
            lss>>read_string;
            if (strcmp(read_string.c_str(), "uniform")== 0){
                solution->soltype = "uniform";
                solution->molins.clear();
                solution->conc.clear();
                solution->concdV.clear();
                int nm;
                lss>>nm;
                for (int i=0; i<nm; i++){
                    double conci;
                    lss>>read_string>>conci;
                    solution->addmol(read_string,conci);
                }
                
                bool comment_found = false;
                while (lss && !comment_found){
                    lss >> read_string;
                    if (strcmp(read_string.c_str(), "Temp")== 0) lss >> solution->Temp;
                    else if (strcmp(read_string.c_str(), "DH_A")== 0) lss >> solution->DH_A;
                    else if (strcmp(read_string.c_str(), "DH_B")== 0) lss >> solution->DH_B;
                    else if (strcmp(read_string.c_str(), "voidV")== 0) lss >> solution->voidV;
                    else if (strcmp(read_string.c_str(), "dV")== 0) lss >> solution->dVtype >> solution->dV;
                    else if (strcmp(read_string.c_str(), "dVvoidV")== 0) lss >> solution->dVvoidV;
                    else if (strcmp(read_string.c_str(), "unitC")== 0) lss >> solution->unitC;
                    else if (strncmp(read_string.c_str(), "#", 1) == 0) {
                        comment_found = true;
                        getout=true;
                    }
                    else {
                        std::string msg = "ERROR: keyword in sol_start  not valid: \""+read_string+"\" \n";
                        error->errsimple(msg);
                    }
                }
                solution->computeNmol();
            }
            else {
                std::string msg = "ERROR: start_sol type not valid. \""+read_string+"\" was found, but only uniform is allowed \n";
                error->errsimple(msg);
            }
            if (me == MASTER) {
                fprintf(screen, "\nDone reading the solution. Temperature is %f and concentrations and numbers for each molecule in box and dV are listed below: \n",solution->Temp);
                for (int i=0; i<chem->Nmol; i++) {
                     fprintf(screen, "  %s   %f  %f %f %f \n",chem->molnames[i].c_str(),chem->mol_cins[i],chem->mol_nins[i],chem->mol_cindV[i],chem->mol_nindV[i]);
                }
            }
        }
        else if (strcmp(word.c_str(), "config0") == 0){
            lss>>inconfig;
            if (strcmp(inconfig.c_str(), "empty") != 0){
                lss >> read_string;
                if (strcmp(read_string.c_str(), "strict") == 0) {
                    time_strict = true;
                    lss >> tstep;
                }
                else if (strcmp(read_string.c_str(), "loose") == 0) {
                    lss >> tstep;
                    time_loose = true;
                }
                else if (strcmp(read_string.c_str(), "first") == 0) {
                    time_first = true;
                }
                else if (strcmp(read_string.c_str(), "last") == 0) {
                    time_last = true;
                }
                else {
                    std::string msg = "ERROR: time type not valid. \""+read_string+"\" was found, but only strict, loose, first, or last, are allowed \n";
                    error->errsimple(msg);
                }
            }
            lss >> read_string;
            if (strcmp(read_string.c_str(), "ru") == 0) {
                lss >> ru;
            }
            else{
                std::string msg = "ERROR: expected argument ru after type of time in config0 command (input or restart file). Insted \""+read_string+"\" was found";
                error->errsimple(msg);
            }
        }
        else if (strcmp(word.c_str(), "fix") == 0){
            std::string read_string3;
            std::getline(lss, read_string3);
            if (me==MASTER) {
                fprintf(screen,"\n inputmsk -- Adding fix:  %s \n",read_string3.c_str());
            }
            fix->add(read_string3);
            
            fprintf(screen,"\n inputmsk -- I am processor %d , part of subcomm %s , and I have %d KMC-free fixes and %d Cont fixes defined now \n",me,(universe->SCnames[universe->color]).c_str(),(int)(fix->fKMCtype.size()),(int)(fix->Ctype.size()));
        }
        else if (strcmp(word.c_str(), "store") == 0){
            std::string read_string3;
            std::getline(lss, read_string3);
            if (me==MASTER) {
                fprintf(screen,"\n inputmsk -- Storing this quantity:  %s \n",read_string3.c_str());
            }
            store->add(read_string3);
            
            fprintf(screen,"\n inputmsk -- Processor %d , part of subcomm %s , has stored: \n",me,(universe->SCnames[universe->color]).c_str());
            for (int i=0; i<store->RegNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->RegCmd[i].c_str());
            for (int i=0; i<store->LatNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->LatCmd[i].c_str());
            for (int i=0; i<store->MinNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->MinCmd[i].c_str());
            for (int i=0; i<store->MinNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->MinModCmd[i].c_str());
        }
        else if (strcmp(word.c_str(), "thermo") == 0){
            if (me==MASTER) msk->wthermo = true;
            lss>>output->th_every;
            lss>>msk->th_fname;
            bool comment_found = false;
            lss >> read_string;
            while (lss && !comment_found){
                if (strncmp(read_string.c_str(), "#", 1) != 0){
                    output->add_thqtt(read_string);
                }
                else  {
                    comment_found = true;
                    getout=true;
                }
                lss >> read_string;
            }
            if (me==MASTER) output->createthermo(msk->th_fname);
        }
        else if (strcmp(word.c_str(), "dump") == 0){
            if (lammpsIO->lammps_active) {
                bool comment_found = false;
                std::string subcname;
                lss >> subcname;
                if (strcmp(subcname.c_str(),(universe->SCnames[universe->color]).c_str())==0 || strcasecmp(subcname.c_str(),"all")==0 ) {
                    std::string read_string,dumpID,dump_group,dump_style,dump_fname,dump_string,tdumpID;
                    int dump_every;
                    lss >> dumpID >> dump_every >> dump_group >> dump_style >> dump_fname;
                
                    // add subcomm name to dump name ----
                    dump_fname = dump_fname+"."+universe->SCnames[universe->color];
                    
                    dump_string = dump_group+" "+dump_style+" 1 "+dump_fname;
                    
                    lss >> read_string;
                    while (lss && !comment_found){
                        if (strncmp(read_string.c_str(), "#", 1) != 0){
                            dump_string = dump_string+" "+read_string;
                        }
                        else  {
                            comment_found = true;
                            getout=true;
                        }
                        lss >> read_string;
                    }
                    output->add_dump(dumpID,dump_every,dump_string);
                }
                else {
                    break;
                }
            }
            else {
                getout=true;
            }
        }
        else if (strcmp(word.c_str(), "tempdump") == 0){
            if (lammpsIO->lammps_active) {
                bool comment_found = false;
                std::string subcname;
                lss >> subcname;
                if (strcmp(subcname.c_str(),(universe->SCnames[universe->color]).c_str())==0 || strcasecmp(subcname.c_str(),"all")==0 ) {
                    std::string read_string, dumpID, dump_string, dump_group, dump_style, dump_every,dump_fname;
                    lss >> dumpID >> dump_group >> dump_style >> dump_every >> dump_fname;
                    dump_fname = dump_fname+"."+universe->SCnames[universe->color];
                    dump_string = "dump "+dumpID+" "+dump_group+" "+dump_style+" "+dump_every+" "+dump_fname;
                    lss >> read_string;
                    while (lss && !comment_found){
                        if (strncmp(read_string.c_str(), "#", 1) != 0){
                            dump_string = dump_string+" "+read_string;
                        }
                        else  {
                            comment_found = true;
                            getout=true;
                        }
                        lss >> read_string;
                    }
                    output->add_tempdump(dumpID,dump_string);
                }
                else {
                    break;
                }
            }
            else {
                getout=true;
            }
        }
        else if (strcmp(word.c_str(), "relax") == 0) {
            if (lammpsIO->lammps_active) {
                bool comment_found = false;
                int every;
                std::string rid,relaxer,rlx_string,rlx_style,rlx_modify,read_string;
                lss >> rid;
                lss >> every;
                lss >> relaxer;
                if (strcmp(relaxer.c_str(), "minimize") == 0){
                    rlx_string = "minimize";
                    for (int i=0;i<4; i++) {
                        lss >> read_string;
                        rlx_string = rlx_string+" "+read_string;
                    }
                    lss >> rlx_style;
                    rlx_modify = "min_modify";
                    lss >> read_string;
                    while (lss && !comment_found){
                        if (strncmp(read_string.c_str(), "#", 1) != 0){
                            rlx_modify = rlx_modify+" "+read_string;
                        }
                        else  {
                            comment_found = true;
                            getout=true;
                        }
                        lss >> read_string;
                    }
                }
                else {
                    std::string msg = "ERROR: only minimize implemented as a relax mode thus far. No nvt, npt, or anything else yet...";
                    error->errsimple(msg);
                }
                
                relax->add_rlx(rid,every,relaxer,rlx_string,rlx_style,rlx_modify);
            }
        }
        else if (strcmp(word.c_str(), "Krun") == 0) {
            double deltat;
            lss >> deltat;
            krun->proceed(deltat);
        }
        else if (strcmp(word.c_str(), "nucleate") == 0) {
            std::string tname, subcom, freqtype;
            freqtype = "step";
            double freqnum = 1;
            double lasteval = 0., starteval=0.;
            lss >> tname >> subcom;
            bool foundend=false, foundevery = false, foundleval=false,foundseval=false;
            while (!foundend) {
                lss >> read_string;
                if (strcmp(read_string.c_str(), "end") == 0 || strcmp(read_string.c_str(), "END") == 0 ) {
                    foundend = true;
                }
                else if (strcmp(read_string.c_str(), "every") == 0){
                    lss >> freqtype >> freqnum;
                    if (strcmp(freqtype.c_str(), "step") != 0 && strcmp(freqtype.c_str(), "time") != 0 ){
                        std::string msg = "ERROR: one nucleation transition has \"every\" option associated to a type of frequency that is neither \"step\" nor \"time\"";
                        error->errsimple(msg);
                    }
                    foundevery = true;
                }
                else if (strcmp(read_string.c_str(), "lasteval") == 0){
                    lss >> lasteval;
                    if (lasteval < 0.) lasteval=0.;
                    foundleval = true;
                }
                else if (strcmp(read_string.c_str(), "starteval") == 0){
                    lss >> starteval;
                    if (starteval < 0.) starteval=0.;
                    foundseval = true;
                }
                
            }
            if (isrestart) {
                if (!foundend) {
                    std::string msg = "ERROR: did not find \"end\" at the end of a nucleation definition in restart file";
                    error->errsimple(msg);
                }
                if (!foundevery) {
                    std::string msg = "ERROR: did not find \"every\" in a nucleation definition in restart file";
                    error->errsimple(msg);
                }
                if (!foundevery) {
                    std::string msg = "ERROR: did not find \"lasteval\" in a nucleation definition in restart file";
                    error->errsimple(msg);
                }
                if (!foundseval) {
                    std::string msg = "ERROR: did not find \"starteval\" in a nucleation definition in restart file";
                    error->errsimple(msg);
                }
                nucleate->addtrans(tname,subcom,freqtype,freqnum,lasteval+freqnum,starteval);
                
            }
            else nucleate->addtrans(tname,subcom,freqtype,freqnum,lasteval,starteval);
            // just a check
            /*
             if (me == MASTER){
             std::string temp_string,temp_string2;;
             temp_string = nucleate->Tnames[nucleate->Ntrans-1];
             temp_string2 = nucleate->scnames[nucleate->Ntrans-1];
             fprintf(screen, "\nNucleation transition \"%s\" added to the preliminary list; later to be associated to subcomm \"%s\"",temp_string.c_str(),temp_string2.c_str());
             }
             */
        }
        else{
            std::string msg = "ERROR: command unknown in input or restart file: "+word;
            error->errsimple(msg);

        }
    }
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Inputmsk::printall()
{
	fprintf(screen,"\n---------ALL ABOUT inputmsk----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
