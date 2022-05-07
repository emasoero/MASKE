#include "spec.h"
#include "solution.h"
#include "chemistry.h"
#include "universe.h"
#include "error.h"
#include <sstream>

using namespace MASKE_NS;

// ---------------------------------------------------------------
// Initialize class
Spec::Spec(MASKE *maske) : Pointers(maske), ipH(-1), ipOH(-1)
{
  IPhreeqc_ptr = new ::IPhreeqc;
}

// ---------------------------------------------------------------
// Class destructor
Spec::~Spec()
{
  delete IPhreeqc_ptr;
}

// ---------------------------------------------------------------
// record new speciation
void Spec::add_spec(std::string id,
                    std::string type,
                    int every,
                    std::string database,
                    const std::vector<std::string>& solvent_names,
                    const std::vector<double>& solvent_molar_masses)
{
  specID.push_back(id);
  spec_type.push_back(type);
  spec_every.push_back(every);
  spec_database.push_back(database);
  if (type == "phreeqc") {
    if (IPhreeqc_ptr->LoadDatabase(database.c_str()) != 0) {
      if (universe->me == 0) IPhreeqc_ptr->OutputErrorString();
    } else {
      if (universe->me == 0) fprintf(screen, "Successfully loaded %s database\n", database.c_str());
    }
  }
  // Find solvents by name
  for (int i = 0; i < solvent_names.size(); i++) {
    int index = -1;
    for (int j = 0; j < chem->Nmol; j++) {
      if (chem->molnames[j] == solvent_names[i])
	index = j;
    }
    if (index > 0) {
      spec_solvents.push_back(index);
    } else {
      error->errsimple("Can't find solvent name " + solvent_names[i] + " in command spec");
    }
  }
  // Copy molar masses
  spec_masses = solvent_molar_masses;
  // Initialise residuals
  for (int i = 0; i < chem->mol_spec.size(); i++) {
    if (chem->mol_spec[i] != "NULL" && chem->mol_spec[i] != "pOH" && chem->mol_spec[i] != "pH") {
      residuals[chem->mol_spec[i]] = 0; // don't actually care if key already exists, just initialise it again
    }
  }
}

// ---------------------------------------------------------------
// write new entry in dump file
void Spec::dospec(int i)
{
  if (spec_type[i] != "phreeqc") {
    if (universe->me == 0) fprintf(screen,"[WARNING] Skipping speciation due to unrecognized type %s\n",spec_type[i].c_str());
    return;
  }
  
  // look for pH and pOH
  if (ipH < 0 || ipOH < 0) {
    ipH = -1;
    ipOH = -1;
    for (int j = 0; j < chem->mol_spec.size(); j++) {
      if (chem->mol_spec[j] == "pOH") {
        if (ipOH > 0 && universe->me == 0) fprintf(screen, "[WARNING] There is more than one molecule flagged as pOH\n");
        ipOH = j;
      } else if (chem->mol_spec[j] == "pH") {
        if (ipH > 0 && universe->me == 0) fprintf(screen, "[WARNING] There is more than one molecule flagged as pH\n");
        ipH = j;
      }
    }
  }
  
  // set up concentration map
  std::map<std::string, double> conc(residuals);
  std::map<std::string, int> mol;    // ids of MASKE molecules associated to phreeqc master species
  for (int j = 0; j < chem->mol_spec.size(); j++) {
    if (chem->mol_spec[j] != "NULL" && chem->mol_spec[j] != "pOH" && chem->mol_spec[j] != "pH") {
      conc[chem->mol_spec[j]] += chem->mol_cins[j];
      mol[chem->molnames[j]] = j;
    }
  }
  // compute mass of water
  double water = 0;
  for (int j = 0; j < spec_solvents.size() ; j++) {
    water += chem->mol_nins[spec_solvents[j]] * spec_masses[j] / nAvo / 1000;
  }
  // set input string for initial concentrations
  std::stringstream ss;
  ss << "SOLUTION 1" << std::endl;
  ss << "units mol/L" << std::endl;
  ss << "pH 7.0 charge" << std::endl; // required to make phreeqc compute pH
  ss << "density 1 calculate" << std::endl;
  ss << "-water " << water << std::endl;
  for (auto it = conc.begin(); it != conc.end(); ++it) {
    ss << it->first << " " << std::scientific << it->second << std::endl; 
  }
  // disable a bunch of outputs all at once
  ss << "SELECTED_OUTPUT" << std::endl;
  ss << "-reset false" << std::endl;
  // set input string for output values
  ss << "USER_PUNCH" << std::endl;
  ss << "10 PUNCH RHO";
  for (auto it = conc.begin(); it != conc.end(); ++it) {
    ss << ", TOT(\"" << it->first << "\")";
  }
  for (auto it = mol.begin(); it != mol.end(); ++it) {
    ss << ", MOL(\"" << it->first << "\")";
  }
  for (auto it = mol.begin(); it != mol.end(); ++it) {
    ss << ", LA(\"" << it->first << "\")";
  }
  ss << ", MOL(\"H+\"), MOL(\"OH-\")" << std::endl;
  std::cout << ss.str() << std::endl;
  // Run phreeqc
  if (this->IPhreeqc_ptr->RunString(ss.str().c_str()) != 0) {
    if (universe->me == 0) IPhreeqc_ptr->OutputErrorString();
  }
  // Get results
  std::vector<VAR> results;
  results.resize(IPhreeqc_ptr->GetSelectedOutputColumnCount());
  for (int j = 0; j < IPhreeqc_ptr->GetSelectedOutputColumnCount(); ++j) {
    VarInit(&results[j]);
    IPhreeqc_ptr->GetSelectedOutputValue(1, j, &results[j]);
  }
  // Transfer results
  int count = 1;
  for (auto it = conc.begin(); it != conc.end(); ++it, ++count) {
    residuals[it->first] = results[count].dVal * results[0].dVal; // TODO: store it in mol and not mol/L
  }
  for (auto it = mol.begin(); it != mol.end(); ++it, ++count) {
    chem->mol_cins[it->second] = results[count].dVal * results[0].dVal;
    residuals[chem->mol_spec[it->second]] -= chem->mol_cins[it->second];
  }
  for (auto it = mol.begin(); it != mol.end(); ++it, ++count) {
    // Get log10 of activities
  }
  // Setting H+ and OH- concentrations
  chem->mol_cins[ipH] = results[count++].dVal * results[0].dVal;
  chem->mol_cins[ipOH] = results[count].dVal * results[0].dVal;
  // Print out debug info
  for (int j = 0; j < chem->mol_spec.size(); j++) {
    std::cout << chem->molnames[j] << ": " << std::scientific << chem->mol_cins[j] << std::endl;
  }
  for (auto it = residuals.begin(); it != residuals.end(); ++it) {
    std::cout << it->first << ": " << std::scientific << it->second << std::endl;
  }  
}

// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Spec::printall()
{
  fprintf(screen,"\n---------ALL ABOUT SPECIATION----------\n");
  fprintf(screen,"---------------------------------------\n\n");
}
