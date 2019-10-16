#include "block.h"
#include <sstream>
#include "universe.h"
#include "error.h"

//#include "lammpsIO.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#include <string.h>

using namespace MASKE_NS;

Block::Block(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    Nreg = 0;
    Nlat = 0;

}



// ---------------------------------------------------------------
// Class destructor
Block::~Block()
{
    
}



// ---------------------------------------------------------------
// Add lattice block to list of existing ones
void Block::add_lattice(std::string instr)
{
    // first check that there is no lattice block already with same name
    for (int i=0; i<LatNames.size(); i++){
        if (strcmp(LatNames[i].c_str(),instr.c_str())==0) {
            std::string msg = "ERROR: lattice block named "+instr+" has been defined more than once\n";
            error->errsimple(msg);
        }
    }
    
    // If OK, add lattice block name and prepare to read through it
    LatNames.push_back(instr);
    tempSvec.clear();
    std::string tstr;
    tstr = "placeholder";
    LatMinCmd.push_back(tstr);
    LatMinMod.push_back(tstr);
  
    MinModFound = false;
    MinCmdFound = false;
    BlockFound = false;
    LatFound = false;
    
    Nlat++;
    
    return;

}



// ---------------------------------------------------------------
// Add region block to list of existing ones
void Block::add_region(std::string instr)
{
    // first check that there is no region block already with same name
    for (int i=0; i<RegNames.size(); i++){
        if (strcmp(RegNames[i].c_str(),instr.c_str())==0) {
            std::string msg = "ERROR: region block named "+instr+" has been defined more than once\n";
            error->errsimple(msg);
        }
    }
    
    // If OK, add lattice region block name and prepare to read through it
    RegNames.push_back(instr);
    tempSvec.clear();
    BlockFound = false;
    RegFound = false;
    
    Nreg++;
    
    return;
}




// ---------------------------------------------------------------
// Add region block to list of existing ones
void Block::add_Latline(std::string instr)
{
    std::istringstream ss(instr);
    std::string word;
    ss >> word;
    
    if (strcmp(word.c_str(), "minimize") == 0){
        LatMinCmd[LatMinCmd.size()-1] = instr;
        MinCmdFound = true;
    }
    else if (strcmp(word.c_str(), "min_modify") == 0){
        LatMinMod[LatMinMod.size()-1] = instr;
        MinModFound = true;
    }
    else if (strcmp(word.c_str(), "block") == 0){
        BlockFound = true;
    }
    else {
        //LatLammps(LatNames.size()).push_back(instr);
        tempSvec.push_back(instr);
        if  (strcmp(word.c_str(), "lattice") == 0){
            LatFound = true;
        }
    }
    
    return;
}



// ---------------------------------------------------------------
// Add region block to list of existing ones
void Block::add_Regline(std::string instr)
{
    //LatLammps(RegNames.size()).push_back(instr);
    std::istringstream ss(instr);
    std::string word;
    ss >> word;
    
    if (strcmp(word.c_str(), "block") == 0){
        BlockFound = true;
    }
    else {
        tempSvec.push_back(instr);
        if  (strcmp(word.c_str(), "region") == 0){
            RegFound = true;
        }
    }
    
    return;
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Block::printall()
{
	fprintf(screen,"\n---------ALL ABOUT BLOCKS----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
