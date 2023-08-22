#include "fix_nufeb.h"
#include "fix.h"
#include "store.h"
#include "lammpsIO.h"
#include "universe.h"
#include "chemistry.h"
#include "error.h"
#include "solution.h"

#include <atom_vec.h>
#include <domain.h>
#include <group.h>
#include <modify.h>

#include <atom_vec_bio.h>
#include <fix_bio_kinetics_ph.h>
#include <fix_bio_kinetics_energy.h>
#include <fix_bio_kinetics_thermo.h>

#include <algorithm>
#include <iomanip>
#include <numeric>

#include <string.h>

using namespace MASKE_NS;

enum Errors { NO_ERROR, KINETICS_NOT_FOUND, DIFFUSION_FOUND };

Fix_nufeb::Fix_nufeb(MASKE *maske) : Pointers(maske)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  init_flag = 0;
  setup_flag = 0;
  setup_exchange_flag = 0;
  kinetics = NULL;
}

// ---------------------------------------------------------------
// Class destructor
Fix_nufeb::~Fix_nufeb() {}

// ---------------------------------------------------------------
// Initialise the fix, see decription in fix_nufeb.h
void Fix_nufeb::init(int pos)
{
  int sid = fix->Csid[pos];

  for (int i=0; i<store->MulCmd[sid].size(); i++) {

    
    lammpsIO->lammpsdo(store->MulCmd[sid][i]);
  }
  init_flag = 1;


  lammpsIO->lammpsdo("timestep 0");
  
  execute(pos, universe->color, 1);
  
}

// ---------------------------------------------------------------
// Setup exchange communication between nufeb and maske 
// subcomm: sub-communicator
void Fix_nufeb::setup(int subcomm)
{
  int error = NO_ERROR;
  if (subcomm == universe->color) {
    LAMMPS_NS::Modify *modify = lammpsIO->lmp->modify;
    for (int i = 0; i < modify->nfix; i++) {
      if (strcmp(modify->fix[i]->style,"kinetics") == 0) {
        kinetics = static_cast<LAMMPS_NS::FixKinetics *>(modify->fix[i]);
      }
    }
    if (!kinetics) error = KINETICS_NOT_FOUND;
  }
  // Checking for error in all processes
  std::vector<int> errors(universe->nprocs);
  MPI_Allgather(&error, 1, MPI_INT, errors.data(), 1, MPI_INT, universe->subcomm);
  for (int i = 0; i < errors.size(); i++) {
    if (errors[i] == KINETICS_NOT_FOUND)
      this->error->errsimple("Fix kinetics not found during fix nufeb setup");
  }
  ncells = kinetics->subn[0] * kinetics->subn[1] * kinetics->subn[2];
  // Check for consistency in NUFEB's molecule parameters for coupling of chemical species concentration.
  // Coupling must strictly be one-to-one with 2 types available: master-to-master, and secondary-to-secondary.
  // Master-to-master coupling associates one master species in MASKE to a master species in NUFEB.
  // This type doesn't account for different chemical forms. Although NUFEB computes each form while
  // computing for speciation only master concentrations will be exchanged.
  // The secondary-to-secondary coupling links two secondary species. MASKE issues a warning if not all forms
  // in NUFEB have an associated species in MASKE. That is to avoid having to track concentration of forms not
  // being considered by MASKE and that must be accounted for when passing on to NUFEB. Moreover, when coupling
  // with phreeqc for speciation, there would be even more complications on how to correctly associate them.
  // To ensure one-to-one associations one might need to create additional chemical species in NUFEB or MASKE,
  // that serve only as buckets for concentrations, or to ensure that charge balance is correct in NUFEB.
  if (universe->key == 0) fprintf(screen,"\nFix NUFEB setup information:\n");
  std::multimap<int, int> mmap;
  for (int i = 0; i < chem->Nmol; i++) {
    int nufeb = chem->mol_nufeb[i];
    // molecule will not be coupled with NUFEB if invalid index range for chemical species is specified
    if (nufeb <= 0 || nufeb > kinetics->bio->nnu + 1) {
      chem->mol_nufeb[i] = -1;
      if (universe->key == 0) fprintf(screen,"Molecule %s is not coupled with NUFEB\n",chem->molnames[i].c_str());
    } else {
      mmap.insert(std::make_pair(nufeb, i));
    }
  }
  const char* forms[] = {"not hydrated", "hydrated", "1st deprotonated", "2nd deprotonated", "3rd deprotonated"};
  for (int i = 1; i <= kinetics->bio->nnu; i++) {
    auto range = mmap.equal_range(i);
    int master = -1;
    // check for master-to-master coupling
    for (auto it = range.first; it != range.second; ++it) {
      int f = chem->mol_nufeb_form[it->second];
      if (f < 0 || f > 4) { // master-to-master coupling
	master = it->second;
	chem->mol_nufeb_form[it->second] = -1;
	if (universe->key == 0) fprintf(screen,"Molecule %s is coupled with NUFEB %s master species\n",chem->molnames[it->second].c_str(),kinetics->bio->nuname[i]);
      }
    }
    // check for secondary-to-secondary coupling
    int form_flags[] = {0, 0, 0, 0, 0};
    for (auto it = range.first; it != range.second; ++it) {
      int f = chem->mol_nufeb_form[it->second];
      if (f >= 0 && f < 5) { // secondary-to-secondary coupling
	if (master >= 0) { // if there is also a master link
	  if (universe->key == 0) fprintf(screen,"[WARNING] Trying to couple with NUFEB %s %s form (%d) in molecule %s is invalid due master coupling in %s molecule\n",kinetics->bio->nuname[i],forms[f],f,chem->molnames[it->second].c_str(),chem->molnames[master].c_str());
	} else {
	  if (kinetics->bio->nugibbs_coeff[i][f] >= INT_MAX) { // pointing to an invalid form
	    if (universe->key == 0) fprintf(screen,"[WARNING] Coupling with invalid NUFEB %s %s form (%d) in molecule %s\n",kinetics->bio->nuname[i],forms[f],f,chem->molnames[it->second].c_str());
	  } else {
	    if (form_flags[f]) { // duplicate link
	      if (universe->key == 0) fprintf(screen,"[WARNING] Duplicate link with NUFEB %s %s form (%d) in molecule %s\n",kinetics->bio->nuname[i],forms[f],f,chem->molnames[it->second].c_str());
	    } else { // finally everything ok
	      if (universe->key == 0) fprintf(screen,"Molecule %s is coupled with NUFEB %s %s form (%d)\n",chem->molnames[it->second].c_str(),kinetics->bio->nuname[i],forms[f],f);
	    }
	  }
	}
      }
      if (kinetics->bio->nugibbs_coeff[i][f] < INT_MAX) form_flags[f] = 1;
    }
    // check if all forms have a molecule link
    if (std::distance(range.first, range.second) > 0 && master < 0) {
      for (int f = 0; f < 5; f++) {
	if (kinetics->bio->nugibbs_coeff[i][f] < INT_MAX && form_flags[f] == 0) {
	  if (universe->key == 0) fprintf(screen,"[WARNING] NUFEB %s %s form (%d) is not associated with a molecule\n",kinetics->bio->nuname[i],forms[f],f);
	}
      }
    }
  }
  // Setup complete
  setup_flag = 1;
}

// ---------------------------------------------------------------
// Setup exchange communication between nufeb and maske 
// subcomm: sub-communicator
void Fix_nufeb::setup_exchange(int subcomm)
{
  buf.resize(1024);
  // exchange subdomain sizes
  sublo.resize(3*universe->nprocs);
  subhi.resize(3*universe->nprocs);
  // Each process will send 3 doubles, pertaining subdomain's lower or upper bounds.
  // sublo vector will contain, in every process, the gathered lower bounds of the subdomains in MPI rank order:
  // (rank0_x, rank0_y, rank0_z, rank1_x, rank1_y, rank1_z, ... , rankn_x, rankn_y, rankn_z)
  MPI_Allgather(lammpsIO->lmp->domain->sublo, 3, MPI_DOUBLE, sublo.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(lammpsIO->lmp->domain->subhi, 3, MPI_DOUBLE, subhi.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);
  // nufeb log output
  if (msk->nulog_flag) {
    msk->nulog << "[rank](domain->sublo): "; 
    for (int i = 0; i < universe->nprocs; i++) {
      msk->nulog << " [" << i << "](" << std::setprecision(9) << sublo[3*i] << "," << sublo[3*i+1] << "," << sublo[3*i+2] << ")";
    }
    msk->nulog << std::endl;
    msk->nulog << "[rank](domain->subhi): "; 
    for (int i = 0; i < universe->nprocs; i++) {
      msk->nulog << " [" << i << "](" << std::setprecision(9) << subhi[3*i] << "," << subhi[3*i+1] << "," << subhi[3*i+2] << ") ";
    }
    msk->nulog << std::endl;
  }
  // check intersections between subdomains
  intersect.resize(universe->nprocs);
  intersect.assign(intersect.size(), true);
  for (int i = 0; i < universe->nprocs; i++) {
    for (int axis = 0; axis < 3; axis++) {
      if (sublo[3*universe->me+axis] >= subhi[3*i+axis] || subhi[3*universe->me+axis] <= sublo[3*i+axis])
        intersect[i] = false;
    }
  }
  // nufeb log output
  if (msk->nulog_flag) {
    msk->nulog << "intersections: "; 
    std::for_each(intersect.begin(), intersect.end(), [&](const int& i) { msk->nulog << " " << i; });
    msk->nulog << std::endl;
  }

  if (subcomm == universe->color) {
    // find out processes to send particles (the ones not belonging to the subcommunicator executing nufeb)
    procs.clear();
    for (int i = 0; i < universe->nprocs; i++) {
      if (universe->color_each[i] != subcomm && intersect[i]) {
        procs.push_back(i);
      }
    }
    // nufeb log output
    if (msk->nulog_flag) {
      msk->nulog << "sending bacteria atoms to ranks: "; 
      std::for_each(procs.begin(), procs.end(), [&](const int& i) { msk->nulog << " " << i; });
      msk->nulog << std::endl;
    }
  } else {
    // find out processes to receive particles from
    procs.clear();
    for (int i = 0; i < universe->nprocs; i++) {
      if (universe->color_each[i] == subcomm && intersect[i]) {
        procs.push_back(i);
      }
    }
    // nufeb log output
    if (msk->nulog_flag) {
      msk->nulog << "receiving bacteria atoms from ranks: "; 
      std::for_each(procs.begin(), procs.end(), [&](const int& i) { msk->nulog << " " << i; });
      msk->nulog << std::endl;
    }
  }
  // Setup complete
  setup_exchange_flag = 1;
}

// ---------------------------------------------------------------
// Compute time increment of the current process
double Fix_nufeb::getDT(int pos)
{
  return fix->Cdt[pos]*fix->Csteps[pos];
}
    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Fix_nufeb::printall()
{
  fprintf(screen,"\n---------ALL ABOUT FIX_NUFEB----------\n");
  fprintf(screen,"---------------------------------------\n\n");
}

// ---------------------------------------------------------------
// Execute nufeb
void Fix_nufeb::execute(int pos, int subcomm, int init)
{
  if (!setup_flag) setup(subcomm);
  // send concentrations to nufeb
  LAMMPS_NS::AtomVecBio *avec = (LAMMPS_NS::AtomVecBio *)lammpsIO->lmp->atom->style_match("bio");
  // vector to store previous negative values of concentrations
  std::vector<double> prev(chem->Nmol, 0);
  // Zero out coupled chemical species concentrations so that we sum them up later on.
  // This accounts for secondary-to-secondary coupling with multiple forms
  for (int i = 0; i < chem->Nmol; i++) {
    prev[i] = std::min(0.0, chem->mol_cins[i]);
    int nufeb = chem->mol_nufeb[i];
    if (nufeb > 0) { // if points to a valid nufeb chemical species
      if (init) {
        for (int j = 0; j < 7; j++) {
          avec->bio->ini_nus[nufeb][j] = 1e-20;
        }
      } else {
        for (int cell = 0; cell < ncells; cell++) {
          kinetics->nus[nufeb][cell] = 1e-20;
        }
      }
    }
  }
  // Sum up contributions of chemical species concentrations
  for (int i = 0; i < chem->Nmol; i++) {
    int nufeb = chem->mol_nufeb[i];
    if (nufeb > 0) { // if points to a valid nufeb chemical species
      if (init) {
	if (chem->mol_cins[i] > 0) {
	  for (int j = 0; j < 7; j++) {
	    avec->bio->ini_nus[nufeb][j] += chem->mol_cins[i];
	  }
        }
      } else {
	if (chem->mol_cins[i] > 0) {
	  for (int cell = 0; cell < ncells; cell++) {
	    kinetics->nus[nufeb][cell] += chem->mol_cins[i];
	  }
        }
      }
    }
  }

  // run nufeb
  std::ostringstream ss1;
  if (!init)
    ss1 << "timestep " << fix->Cdt[pos];
  lammpsIO->lammpsdo(ss1.str());
  std::ostringstream ss2;
  ss2 << "run " << fix->Csteps[pos];
  if (!init)
    ss2 << " pre no";
  ss2 << " post no";
  lammpsIO->lammpsdo(ss2.str());
  lammpsIO->lammpsdo("timestep 0");
  // Compute the difference in concentrations after running NUFEB
  std::vector<double> dconc(chem->Nmol, 0);
  std::vector<double> temp_buc(chem->Nmol, 0); //a temporary bucket to store the sum of concentrations from other forms
  for (int i = 0; i < chem->Nmol; i++) {
    int nufeb = chem->mol_nufeb[i];
    if (nufeb > 0) { // if points to a valid nufeb chemical species
      double conc = 0;
      if (universe->key == 0) {
	      int form = chem->mol_nufeb_form[i];
	      if (form >= 0) { // secondary-to-secondary coupling
	        dconc[i] = kinetics->activity[nufeb][form][0] + prev[i] - chem->mol_cins[i];

          for (int spc=0;spc<5;spc++) fprintf(screen,"Concentration of form (%d) is %f \n",spc,kinetics->activity[nufeb][spc][0]);
          sleep(1);

          fprintf(screen,"Total concentration of the nutrient (%d) is %f \n",nufeb,kinetics->nus[nufeb][0]);
          sleep(1);
    
          temp_buc[i]=kinetics->nus[nufeb][0]-kinetics->activity[nufeb][form][0];
          fprintf(screen,"The difference in concentration of the nutrient (%d) and the form (%d) is %f \n",nufeb,form, temp_buc[i]);
          sleep(1);


	      } 
        else { // master-to-master coupling
	          dconc[i] = kinetics->nus[nufeb][0] + prev[i] - chem->mol_cins[i];
	      }
      }
    }
  }
  // Compute new concentrations based on the difference
  solution->updateconc(pos, dconc);
  // Send the new total concentrations from nufeb
  for (int i = 0; i < chem->Nmol; i++) {
    int nufeb = chem->mol_nufeb[i];
    if (nufeb > 0) { // if points to a valid nufeb chemical species
      //fprintf(screen,"\n\n DEBUG 1: Proc %d, concentration of nutrient (%d) mapped to molecule (%d) is %e \n",me, nufeb, i, chem->mol_cins[i]);
      //sleep(1);

      if (!init){
         if (universe->key == 0){
            for (int j=0; j<universe->nprocs;j++){
            if (j!=me) {
              MPI_Send(&chem->mol_cins[i],1, MPI_DOUBLE,j,j,MPI_COMM_WORLD);
              MPI_Send(&chem->mol_nins[i],1, MPI_DOUBLE,j,j,MPI_COMM_WORLD);
              MPI_Send(&chem->mol_cindV[i],1, MPI_DOUBLE,j,j,MPI_COMM_WORLD);
              MPI_Send(&chem->mol_nindV[i],1, MPI_DOUBLE,j,j,MPI_COMM_WORLD);
            }
            }
          }
          else {
            MPI_Recv(&chem->mol_cins[i],1, MPI_DOUBLE,universe->subMS[subcomm],me,MPI_COMM_WORLD,&status);
            MPI_Recv(&chem->mol_nins[i],1, MPI_DOUBLE,universe->subMS[subcomm],me,MPI_COMM_WORLD,&status);
            MPI_Recv(&chem->mol_cindV[i],1, MPI_DOUBLE,universe->subMS[subcomm],me,MPI_COMM_WORLD,&status);
            MPI_Recv(&chem->mol_nindV[i],1, MPI_DOUBLE,universe->subMS[subcomm],me,MPI_COMM_WORLD,&status);
          }
      }
     
      //if (!init) MPI_Bcast(&chem->mol_cins[i], 1, MPI_DOUBLE, universe->subMS[subcomm], MPI_COMM_WORLD);

    }
  }
}

// ---------------------------------------------------------------
// Exchange bacteria atoms to maske subcommunicators
// id: nufeb continuous process id
// subcomm: sub-communicator
void Fix_nufeb::exchange(int id, int subcomm)
{
  if (!setup_exchange_flag) setup_exchange(subcomm);
  if (subcomm == universe->color) {
    // pack atoms inside the subdomain of other processes
    nsend.resize(procs.size());
    nsend.assign(nsend.size(), 0);
    int group = lammpsIO->lmp->group->find(fix->Cgroups[id].c_str());
    int bitmask = lammpsIO->lmp->group->bitmask[group];
    LAMMPS_NS::Atom *atom = lammpsIO->lmp->atom;
    int *mask = atom->mask;
    double **x = atom->x;
    std::vector<MPI_Request> requests(procs.size());
    int total = 0; // total number of doubles to be sent
    for (int p = 0; p < procs.size(); p++) {
      int r = procs[p]; // process rank
      for (int i = 0; i < atom->nlocal; i++) {
        if (x[i][0] >= sublo[3*r] && x[i][0] < subhi[3*r]
            && x[i][1] >= sublo[3*r+1] && x[i][1] < subhi[3*r+1]
            && x[i][2] >= sublo[3*r+2] && x[i][2] < subhi[3*r+2]) {
          if (mask[i] & bitmask) {
            if (buf.size() < total + 1024) {
              buf.resize(2*buf.size());
            }
            int n = atom->avec->pack_exchange(i,&buf[total]);
            nsend[p] += n;
            total += n;
          }
        }
      }
      // send the number of doubles to be sent
      MPI_Isend(&nsend[p], 1, MPI_INT, r, 0, MPI_COMM_WORLD, &requests[p]);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    // nufeb log output
    if (msk->nulog_flag) {
      msk->nulog << "[rank](number of doubles to send): ";
      for (int p = 0; p < procs.size(); p++) {
        msk->nulog << "[" << p << "](" << nsend[p] << ") ";
      }
      msk->nulog << std::endl;
    }
    // exchange atoms
    int begin = 0;
    for (int p = 0; p < procs.size(); p++) {
      int r = procs[p]; // process rank
      MPI_Send(&buf[begin], nsend[p], MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
      begin += nsend[p];
    }
  } else if (lammpsIO->lammps_active) {
    // Receive the number of doubles to be sent
    nrecv.resize(procs.size());
    nrecv.assign(nrecv.size(), 0);
    std::vector<MPI_Request> requests(procs.size());
    for (int i = 0; i < procs.size(); i++) {
      MPI_Irecv(&nrecv[i], 1, MPI_INT, procs[i], 0, MPI_COMM_WORLD, &requests[i]);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    // nufeb log output
    if (msk->nulog_flag) {
      msk->nulog << "[rank](number of doubles to receive): "; 
      for (int i = 0; i < procs.size(); i++) {
        msk->nulog << "[" << procs[i] << "](" << nrecv[i] << ") ";
      }
      msk->nulog << std::endl;
    }
    // ensure buf is big enough to host incoming data
    int total = std::accumulate(nrecv.begin(), nrecv.end(), 0); // total number of doubles to receive
    if (buf.size() < total)
      buf.resize(1.5*total);
    // exchange atoms
    int begin = 0;
    for (int p = 0; p < procs.size(); p++) {
      int r = procs[p]; // process rank
      MPI_Irecv(&buf[begin], nrecv[p], MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &requests[p]);
      begin += nrecv[p];
    }
    if (requests.size() > 0) MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    LAMMPS_NS::Atom *atom = lammpsIO->lmp->atom;
    // delete bacteria atoms
    std::stringstream ss;
    ss << "delete_atoms group " << fix->aCgroups[id] << " compress no";
    lammpsIO->lammpsdo(ss.str());
    // unpack received atoms
    int m = 0;
    while (m < total) {
      m += atom->avec->unpack_exchange(&buf[m]);
      atom->tag[atom->nlocal-1] = 0;
    }
    // update atom->natoms: total number of atoms
    LAMMPS_NS::bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,universe->subcomm);
    // update map
    atom->tag_extend();
    if (atom->map_style) {
      atom->map_init(1);
      atom->map_set();
    }
  }
}
