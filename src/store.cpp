#include "store.h"
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

Store::Store(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    Nreg = 0;
    Nlat = 0;
    Nmin = 0;
}



// ---------------------------------------------------------------
// Class destructor
Store::~Store()
{
    
}



// ---------------------------------------------------------------
// Add lattice block to list of existing ones
void Store::add(std::string instr)
{
    
    std::istringstream ss(instr);
    std::string word;
    ss >> word;
    
    if (strcmp(word.c_str(), "region") == 0){
        ss >> word;
        RegNames.push_back(word);
        RegCmd.push_back(instr);
        Nreg++;
        //fprintf(screen,"Region %s command stored: %s\n",RegNames[Nreg-1].c_str(),RegCmd[Nreg-1].c_str());
    }
    else if (strcmp(word.c_str(), "lattice") == 0){
        ss >> word;
        LatNames.push_back(word);
        std::string read_string;
        std::getline(ss, read_string);
        read_string = "lattice "+read_string;
        LatCmd.push_back(read_string);
        LatDVnames.push_back("none");
        LatDVcmd.push_back("none");
        Nlat++;
        //fprintf(screen,"Lattice %s command stored: %s\n",LatNames[Nlat-1].c_str(),LatCmd[Nlat-1].c_str());
    }
    else if (strcmp(word.c_str(), "DV") == 0){
        ss >> word;
        bool found = false;
        for (int i=0; i<LatNames.size(); i++) {
            if (strcmp(word.c_str(), LatNames[i].c_str())==0) {
                found = true;
                ss >> word;   //the word variable: not to be used actually but just as a memo for the user..
                ss >> LatDVnames[i];   // the name of DV
                std::string read_string;
                std::getline(ss, read_string);
                LatDVcmd[i] = read_string;
            }
        }
        if (found==false) {
            std::string msg = "ERROR: a DV could not be stored because lattice "+word+" was not found\n";
            error->errsimple(msg);
        }
    }
    else if (strcmp(word.c_str(), "minimize") == 0){
        ss >> word;
        MinNames.push_back(word);
        std::string read_string;
        read_string = "minimize";
        for (int i=0; i<4; i++) {
            ss >> word;
            read_string = read_string + " " + word;
        }
        MinCmd.push_back(read_string);
        std::getline(ss, read_string);
        MinModCmd.push_back(read_string);
        Nmin++;
        //fprintf(screen,"Minimiser %s command stored: %s with modifier %s \n",MinNames[Nmin-1].c_str(),MinCmd[Nmin-1].c_str(),MinModCmd[Nmin-1].c_str());
    }
    else {
        std::string msg = "ERROR: store can only take region, lattice, or minimize as input. Instead   "+word+" was found\n";
        error->errsimple(msg);
    }
    
    return;

}




/*

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

*/


// ---------------------------------------------------------------
// Printing info about the Store class (possibly useful for debugging)
void Store::printall()
{
	fprintf(screen,"\n---------ALL ABOUT STORE----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
