#include "fix_nufeb.h"
#include "fix.h"
#include "store.h"
#include "lammpsIO.h"
#include "universe.h"
#include "chemistry.h"
#include "error.h"

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

using namespace MASKE_NS;

enum Errors { NO_ERROR, KINETICS_NOT_FOUND, DIFFUSION_FOUND };

Fix_nufeb::Fix_nufeb(MASKE *maske) : Pointers(maske)
{
  setup_flag = 0;
  setup_exchange_flag = 0;
  first_flag = 1;
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
  ncells = kinetics->nx * kinetics->ny * kinetics->nz;
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
  // (rank0_x, rank0_y, rank0_z, rank1_x, rank2_y, rank1_z, ... , rankn_x, rankn_y, rankn_z)
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
void Fix_nufeb::execute(int pos, int subcomm)
{
  if (!setup_flag) setup(subcomm);
  // send concentrations to nufeb
  LAMMPS_NS::AtomVecBio *avec = (LAMMPS_NS::AtomVecBio *)lammpsIO->lmp->atom->style_match("bio");
  for (int i = 0; i < chem->Nmol; i++) {
    int nufeb = chem->mol_nufeb[i];
    if (nufeb > 0) { // if points to a valid nufeb chemical species
      if (first_flag) {
	for (int j = 0; j < 7; j++) {
	  avec->bio->ini_nus[nufeb][j] = (chem->mol_cins[i] < 0) ? 1e-20 : chem->mol_cins[i];
	}
      } else {
	for (int cell = 0; cell < ncells; cell++) {
	  kinetics->nus[nufeb][cell] = (chem->mol_cins[i] < 0) ? 1e-20 : chem->mol_cins[i];
	}
      }
    }
  }

  // run nufeb
  std::ostringstream ss1;
  ss1 << "timestep " << fix->Cdt[pos];
  lammpsIO->lammpsdo(ss1.str());
  std::ostringstream ss2;
  ss2 << "run " << fix->Csteps[pos] << " pre no post no";
  lammpsIO->lammpsdo(ss2.str());
  // Send the new total concentrations from nufeb 
  for (int i = 0; i < chem->Nmol; i++) {
    int nufeb = chem->mol_nufeb[i];
    if (nufeb > 0) { // if points to a valid nufeb chemical species
      double conc;
      if (universe->key == 0) {
	conc = kinetics->nus[nufeb][0];
      }
      MPI_Bcast(&conc, 1, MPI_DOUBLE, universe->subMS[subcomm], MPI_COMM_WORLD);
      chem->mol_cins[i] = conc;
    }
  }
  first_flag = 0;
}

// ---------------------------------------------------------------
// Execute nufeb
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
  } else {
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
