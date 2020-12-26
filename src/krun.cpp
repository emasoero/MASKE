#include "krun.h"
#include "fix.h"
#include "randm.h"
#include "fix_delete.h"
#include "fix_nucleate.h"
#include "fix_Cfoo.h"
#include "universe.h"
#include "lammpsIO.h"
#include "output.h"
#include "solution.h"
#include "chemistry.h"
#include "relax.h"
#include "error.h"

#ifdef MASKE_WITH_NUFEB
#include "fix_nufeb.h"
#endif
#ifdef MASKE_WITH_SPECIATION
#include "spec.h"
#endif

#include <algorithm>
#include <numeric>
#include <iomanip>

#include <group.h>
#include <atom_vec.h>
#include <memory.h>
#include <domain.h>
#include <comm.h>
#include <modify.h>
#include <fix.h>

#include <string.h>
 
using namespace MASKE_NS;

union ubuf {
  double d;
  int64_t i;
  ubuf(double arg) : d(arg) {}
  ubuf(int64_t arg) : i(arg) {}
  ubuf(int arg) : i(arg) {}
};

Krun::Krun(MASKE *maske) : Pointers(maske)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    reset_event_vec = true;   // if the simulation is not a restart, then at the first step the code will compute the number of events in each RF-KMC process and prepare empty vectors to accommodate rates coming from the fixes
    //Ntrans = 0;
    //Ntsc=0;
    QkTint = 0.;
    resetQkTint = true;

    nufeb_buf.resize(1024);
}



// ---------------------------------------------------------------
// Class destructor
Krun::~Krun()
{
    
}



// ---------------------------------------------------------------
// Class destructor
void Krun::proceed(double deltat)
{
    
    // screen output to debug
     int sleeptime = me;    // KMR - defines an integer to give as input to sleep() function just used for debigging.   Gives numreical valuw to integer equal to rank of current processor in COMM_WORLD
    /*
     sleep(sleeptime);
    fprintf(screen,"\n\n Proc %d, SUBCOM %d Ready to start krun \n",me,universe->color);
    MPI_Barrier(MPI_COMM_WORLD);
    if (me==MASTER) {
        fprintf(screen,"\n\n Let's start krun \n");
    }
    sleep(2);
     */
    
    // some initial variables and quantities
    double start_time = msk->tempo;   // the time at which this krun loop started
    double end_time = start_time + deltat;   //time past/reached which this loop stops
    double u;       // Random number  for later usage, to compute  future occurrence time (FOT) of  RF-KMC events
    if (me==MASTER) u = randm->pick01();
    
    
    
    if (resetQkTint){
        QkTint = 0.;
        resetQkTint = false; // by default, the next krun will inherit leftover cumulated rate from current run, unless something is defined in between runs which requires QKTint to be reset
    }
    bool write_first_thermo = true;
    bool write_first_dumps = true;

    // Setting time step to 0
    lammpsIO->lammpsdo("timestep 0");

    for (int i=0; i<fix->Ctype.size(); i++) {
#ifdef MASKE_WITH_NUFEB
      if (fix->Ctype[i] == "nufeb") {
	fix_nufeb->init(i);
      }
#endif
    }
    
    //**********************************
    // FROM NOW ON, TIME WILL ADVANCE
    while (msk->tempo < end_time){
        
        
        // just to mix up lammps config a bit and force lammps to use multiple processors..
        
         /*
          if (lammpsIO->lammps_active){
            std::string tolmp = "displace_atoms all random 0.0050 0.0050 0.0080 123224";
            //std::string tolmp = "displace_atoms all move 10 20 -4";
            lammpsIO->lammpsdo(tolmp);
            tolmp = "run 1";   // just to see something in dump... to be removed later
            lammpsIO->lammpsdo(tolmp);
        }
          */
        
        
        
        if (write_first_thermo) {
	    output->writethermo();
            write_first_thermo = false;
        }
        
        if (write_first_dumps) {
            for (int i = 0; i<output->dumpID.size();i++){
	        output->writedump(i);
                output->dump_first[i]=false;  // from now on, not the first time the dump is written anymore
            }
            write_first_dumps = false;
        }
        
        
        if (me==MASTER) fprintf(screen,"\n MASKE tempo is %e \n  start_time was %f \n End time is %f",msk->tempo,start_time,end_time);
        
        // screen output to debug
        //sleep(me);
        fprintf(screen,"\n\n Proc %d, SUBCOM %d entering step at tempo %e \n",me,universe->color,msk->tempo);
        MPI_Barrier(MPI_COMM_WORLD);
        if (me==MASTER) fprintf(screen,"\n\n Let's start the tempo %e iteration \n",msk->tempo );
       // sleep(1);
    
        
        
        //-----------------------------------------
        // CLEAR ALL EVENT VECTORS AND ASSIGN SIZES TO NEW ONES
        if (reset_event_vec) {
            
            QkTint = 0.;
            
            // screen output to debug
            //sleep(me);
            fprintf(screen,"\n\n Proc %d, SUBCOM %d tempo %e, resetting vectors \n",me,universe->color,msk->tempo);
            MPI_Barrier(MPI_COMM_WORLD);
            //sleep(1);
            
            fix_del->clearvecs();
            fix_nucl->clearvecs();   // You can consider moving here the deletion of all existing trial particles for nucleation - now it is done fix by fix in the init function below (which works too)
            
            
            // each processor scans through the list of their RF-KMC processes
            int tot_DEL_evt = 0;
            int tot_NUCL_evt = 0;
            for (int i=0; i<fix->fKMCtype.size(); i++) {
        
                if (strcmp(fix->fKMCtype[i].c_str(),"delete")==0) {
                    fix_del->init(i);   // counts number of particles deletable by current fix, and records it in position i of fix->fKMCnevents vector
                    fix->fKMCfirst[i] = tot_DEL_evt;
                    int Ng = fix->fKMCnevents[i];
                    tot_DEL_evt = tot_DEL_evt + Ng;
                    fix_del->add_Ng_apos(Ng,fix->fKMCaID[i]);  // add number of deletable particles to vector cumulating it for all fixes, and absolute position of fix in corresponding vector
                }
                else if (strcmp(fix->fKMCtype[i].c_str(),"nucleate")==0) {
                    fix_nucl->delete_trial(i);  // delete trial particles in this fix, if any exists
                    fix_nucl->init(i); // creates trial particles for this fix, count them, and record the count inin position i of fix->fKMCnevents
                    fix->fKMCfirst[i] = tot_NUCL_evt;
                    int Ng = fix->fKMCnevents[i];
                    tot_NUCL_evt = tot_NUCL_evt + Ng;
                    fix_nucl->add_Ng_apos(Ng,fix->fKMCaID[i]);  // add number of trial particles to vector cumulating it for all fixes, and absolute position of fix in corresponding vector
                }
                // same to be added for any other KMC event...
            }

            fMC_justreset = true;   // to tell the fix whether new rate vectors need be created or existing ones to be read
        }
        reset_event_vec = false;   // to avoid resetting vectors again and next step, unless something happens to turn the reset on again
        //-----------------------------------------

        
        
        
        
        //-----------------------------------------
        // SAMPLE RF-KMC EVENTS AND COMPUTE THEIR FUTURE OCCURRENCE TIME (FOT)
        bool rfKMCcrit = true;  //just a placeholder for now. In the future the criteria will be in a vector in fix.cpp, associated to the ith process
        
        if (rfKMCcrit) {
            // each processor scans through the list of their RF-KMC processes
            for (int i=0; i<fix->fKMCtype.size(); i++) {
                if (strcmp(fix->fKMCtype[i].c_str(),"delete")==0)           fix_del->sample(i);
                else if (strcmp(fix->fKMCtype[i].c_str(),"nucleate")==0)    fix_nucl->sample(i);
            }
        }
        
        double Qk = 0.;   // cumulative rate
        double Qkdel = 0.;   // cumulative rate from delete events only
        double Qknucl= 0.;   // cumulative rate from nucleation events only
        
        if (universe->key==0) {
            Qkdel = fix_del->CrateAll();   //extract total rate of all delete and events in subcomm
            Qknucl = fix_nucl->CrateAll();
        }
        Qk = Qkdel + Qknucl;  //  + ... all RF-KMC event types to be defined in future
        
        
        // each submaster sends its Qk to the master (me==0)
        double *Qkall;  // an array for the master to record the Qk from each submaster
        Qkall = new double[universe->nprocs];
        for (int i=0; i<universe->nprocs; i++) Qkall[i] = 0.;
        
        if (universe->key==0) {
            Qkall[me] = Qk;  //each submaster records its Qk into its local Qkall, which is now ready to be passed to the global MASTER
            fprintf(screen,"\n\n Vector of cumulative delete + nucleate rates from processor %d is: \n",me);
            for (int i=0; i<universe->nprocs; i++) {
                fprintf(screen," %e ",Qkall[i]);
            }
        }

        if (me>0){
            int dest = 0;
            MPI_Send(&Qkall[me], 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
        }
        if (me == MASTER) {
            for (int source=1; source<universe->nprocs; source++) {
                MPI_Recv(&Qkall[source], 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
            }
            fprintf(screen,"\n\n Vector of cumulative delete + nucleate rates from all subcomms is: \n");
            for (int i=0; i<universe->nprocs; i++) {
                fprintf(screen," %e ",Qkall[i]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Master sums up Qks received from all submasters into the total Qktot at this time
        double Qktot = 0.;
        if (me == MASTER) {
            for (int i=0; i<universe->nprocs; i++) {
                Qktot += Qkall[i];
            }
        }
        
        delete [] Qkall;
        
        // Master computes waiting time of next RF-KMC event
        double DtKMC = 0.;
        if (me == MASTER) {
            // KMC condition is rate x Dt = log(1/u)
            // here the rate x Dt is split between already cumulated one (QkTint) and future one, Qktot x DtKMC (the latter being the time from now until the next event, not the time from previous to next event: some of the time from the previous event has already been accumulated in QkTint)
            // So the generating equation is QkTint + DtKMC x Qktot = log(1/u), which is solved for DtKMC
            DtKMC = ( log(1./u) - QkTint ) / Qktot;
            fprintf(screen,"\n\n DtKMc at tempo %e is: %e\n",msk->tempo,DtKMC);
        }
        
        
        
            
        // Compute waiting time for all continuus processes
        
        int *nCpSC;   // array storing the number of continuum processes in each procesor (i.e. in subcomm to which processor pertains)
        nCpSC = new int [universe->nprocs];
        for (int i=0; i<universe->nprocs; i++) nCpSC[i] = 0;
        nCpSC[me] = fix->Ctype.size();
        
        int *pcol;   // array with color of subcomm of each proc (-1 for all except current proc)
        pcol = new int [universe->nprocs];
        for (int i=0; i<universe->nprocs; i++) pcol[i] = -1;
        pcol[me] = universe->color;
        
        
        double *DtC; // array storing the Dts of all continuum processes in the current subcom. The array starts always with a -1, to avoid communication problems for subcomm with 0 cont processes defined for them (-1 if proc is slave, >=0 for submasters)
        DtC = new double [fix->Ctype.size()+1];  // the +1 is to accommodate the initial -1 dummy entry
        for (int i=0; i<fix->Ctype.size()+1; i++) DtC[i]=-1;
        
        bool Ccrit = true;  //just a placeholder for now. In the future the criteria will be in a vector in fix.cpp, associated to the ith process
        if (Ccrit) {
            // each processor scans through the list of their Cont processes
            for (int i=0; i<fix->Ctype.size(); i++) {
                if (strcmp(fix->Ctype[i].c_str(),"foo")==0) {
                    double temp_Dt = fix_cfoo->getDT(i);     // +1 because I want DtCall[0]=-1. If static DT just read it from vector in fix.cpp. If dynamic DT, it may need to be computed somhow, hence I defined a placeholder function.  Each proc in subcom knows about Dt, but we will only work with the value in the submaster..
                    if (universe->key==0) {
                        DtC[i+1] = fix->Cleval[i] + temp_Dt - msk->tempo;
                    }
                }
#ifdef MASKE_WITH_NUFEB
		if (strcmp(fix->Ctype[i].c_str(),"nufeb")==0) {
		  double temp_Dt = fix_nufeb->getDT(i);
		  if (universe->key==0) {
		    DtC[i+1] = fix->Cleval[i] + temp_Dt - msk->tempo;
		  }
		}
#endif
            }
        }
        
        // Master gatehrs number of Cont processes for each processor
        if (me > 0) {
            int dest = 0;
            MPI_Send(&nCpSC[me], 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&pcol[me], 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
           
        }
        if (me==MASTER){
            for (int source=1; source<universe->nprocs; source++) {
                MPI_Recv(&nCpSC[source], 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&pcol[source], 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
            }
            fprintf(screen,"\n\n Number of continuous proceeses gathered by master \n Proc_No   Subcom   nCprocs\n");
            for (int i=0; i<universe->nprocs; i++) {
                fprintf(screen,"%d %s %d\n",i,universe->SCnames[pcol[i]].c_str(),nCpSC[i]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Master gatehrs Dt arrays for each Cont processes from each processor.   The size of each array is 1+nCpSC[me]
        double * DtCall;   // array in master assemblying the processor-specific DtC arrays
        int *CplocID;   // array in master containing local id of each cont process in DtCall in its own subcom
        int *pcolAll;   // array in master with color (subcomm id) of process in each entry in DtCall
        int Dtsize = 1;
        if (me==MASTER) {
            Dtsize = 0;
            for (int i=0; i<universe->nprocs; i++) Dtsize += 1 +nCpSC[i];
            DtCall = new double [Dtsize];
            CplocID = new int [Dtsize];
            
            // writing color (subcom ID) of all processes whose Dt is recorded in  DtCall
            pcolAll = new int [Dtsize];
            int ppos = 0;
            for (int i=0; i<universe->nprocs; i++) {
                for (int j=0; j< 1+nCpSC[i]; j++) {
                    pcolAll[ppos] = pcol[i];
                    ppos++;
                }
            }
        }
       
        
        // master gathers DtC arrays from each processor (meaningful only from submasters) inside array DtCall.
        // While doing this, master also writes local id of each process in DtCall, filling array CplocID
        if (me > 0) {
            int dest = 0;
            MPI_Send(&DtC[0], 1+nCpSC[me], MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            
        }
        if (me==MASTER){
            for (int i=0; i<nCpSC[me]+1; i++) {
                DtCall[i] = DtC[i];
            }
            int ppos = 0;    //position where to record DtCall array
            for (int source=1; source<universe->nprocs; source++) {
                ppos += nCpSC[source-1]+1;
                MPI_Recv(&DtCall[ppos], 1+nCpSC[source], MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
            }
            ppos = 0;
            fprintf(screen,"\n\n Dts of all continuous proceeses gathered by master \n Proc_No   Subcom   DtC   proc_ID_loc\n");
            for (int i=0; i<universe->nprocs; i++) {
                for (int j=0; j< 1+nCpSC[i]; j++) {
                    CplocID[ppos]=j-1;
                    fprintf(screen,"%d %s %f %d\n",i,universe->SCnames[pcolAll[ppos]].c_str(),DtCall[ppos],CplocID[ppos]);
                    ppos++;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        delete [] DtC;
        delete [] nCpSC;
        delete [] pcol;
        
        
        
        
        // Master selects smallest FOT among all KMC and continuous types above  --> next event type
        double Dt = 0.;
        int KMCexecute = 0;  // a bool actually, 0 = false, 1 = true
	int Cexecute = 0;    // flags the execution of a countinous process
        int Cpid2exec = -1;  // local ID of continuous process to execute
        int SC2exec = -1;    // ID of subomm supposed to carry out the continuous event
        int endofrun = 0;    // a bool actually, 0 = false, 1 = true
	std::string type;
	
        if (me == MASTER) {
            Dt = DtKMC;
            if (Qktot > 0.) {
	      KMCexecute = 1;
	      Cexecute = 0;
	    }
            if (Dt < 0.)  Dt = 0.;    // the KMC process is chosen because associated to time < tempo.. acceptable.. but a warning should be generated..
            else {
                for (int i=0; i<Dtsize; i++) {
                    if (DtCall[i] >= 0. && DtCall[i]<Dt) {
                        Dt = DtCall[i];  //the corresponding subcom color is in position i of pcolAll array
                        KMCexecute = 0;
			Cexecute = 1;
                        Cpid2exec = CplocID[i];
                        SC2exec = pcolAll[i];
			type = fix->aCtype[Cpid2exec];
                    }
                }
            }
            if (Dt > end_time - msk->tempo){
                Dt = end_time - msk->tempo;   // This means that loop is over: all continuum processes should be executed in this case and QkTint of KMC events should be updated
                endofrun = 1;
            }
        }
        
        if (me==MASTER) {
            delete [] DtCall;
            delete [] CplocID;
            delete [] pcolAll;
        }
        
        
        // Master communicates selected dt and type of chosen process to all processors
        if (me==MASTER){
            for (int dest=1; dest<universe->nprocs; dest++) {
                MPI_Send(&Dt, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
                MPI_Send(&KMCexecute, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
                MPI_Send(&Cexecute, 1, MPI_INT, dest, 3, MPI_COMM_WORLD);
                MPI_Send(&Cpid2exec, 1, MPI_INT, dest, 4, MPI_COMM_WORLD);
                MPI_Send(&SC2exec, 1, MPI_INT, dest, 5, MPI_COMM_WORLD);
                MPI_Send(&endofrun, 1, MPI_INT, dest, 6, MPI_COMM_WORLD);
		MPI_Send(type.c_str(), type.length()+1, MPI_CHAR, dest, 7, MPI_COMM_WORLD);
            }
        }
        if (me > 0) {
            int source = 0;
            MPI_Recv(&Dt, 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&KMCexecute, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&Cexecute, 1, MPI_INT, source, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(&Cpid2exec, 1, MPI_INT, source, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(&SC2exec, 1, MPI_INT, source, 5, MPI_COMM_WORLD, &status);
            MPI_Recv(&endofrun, 1, MPI_INT, source, 6, MPI_COMM_WORLD, &status);
	    MPI_Status status;
	    MPI_Probe(source, 7, MPI_COMM_WORLD, &status);
	    int length;
	    MPI_Get_count(&status, MPI_CHAR, &length);
	    char *type_recv = new char[length];
	    MPI_Recv(type_recv, length, MPI_CHAR, source, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    type = std::string(type_recv);
	    delete [] type_recv;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        // all processors update tempo
        msk->tempo = msk->tempo + Dt;
        if (me == MASTER) {
            fprintf(screen,"\n\n New tempo is %e \n",msk->tempo);
        }
        
        
        // all rates of KMC processes are integrated by Dt by all processors. Master integrates the grand cumulate QkTint too
        double *QkT, *QkTC;  // position 0 refers to fix_del, position 1 to fix_nucleate.
        int nKMCs = 2;      // number of KMC fix types allowed in the program. For now, delete and nucleate, so 2.
        QkT = new double[nKMCs];    // cumulated time integrated rates of all fix_KMCs in this  proc (nonzero only if submaster)
        QkTC = new double[nKMCs];    // sum of QkT above over fix_del and fix_nucl in same subcomm
        // NB: if you have subcomm A and B, with fix_delA1 and fix_delA2 in A and fix_delB1 and fix_nuclB1 in B, QkT for A will be QkT[0] = sum of fix_delA1+A2, QkT[1] = 0 (cause no fix_nucl in A), whereas for B it will be QkT[0] = fix_delB1 and QkT[1] = fix_nuclB1.      Then, QkTC be QkTC[0] = QkT[0] and QkTC[1] = QkT[0] + Qkt[1] for either of the subcomms.
        
        //if (msk->wplog) { std::string msg = "\n Eddig jo\n"; output->toplog(msg);}
        
        double *QkTCall, *QkTCallC;    // arrays to be used by master to gather QkTC[end] values from all processors, and to cumulate them
        if (me==MASTER) {
            QkTCall = new double[universe->nprocs];
            QkTCallC = new double[universe->nprocs];
        }
        
        for (int i=0; i<nKMCs; i++) QkT[i] = 0.;
        if (universe->key==0) {
             QkT[0] = fix_del->CratesTupdate(Dt);   //updates the vector of time integrals of cumulated rates, CratesT, in fix_delete and fix_nucleate
             QkT[1] = fix_nucl->CratesTupdate(Dt);
        }
        QkTC[0] = QkT[0];
        for (int i=1; i<nKMCs; i++) QkTC[i] = QkTC[i-1] + QkT[i];

        if (msk->wplog) {
            std::string msg = "\n All-fix cumrate vectors to send to master: \n QkT \t QkTC\n";
            std::ostringstream ss;
            for (int i=0; i<nKMCs; i++) {
                ss << QkT[i];   msg += ss.str()+"\t";  ss.str("");   ss.clear();
                ss << QkTC[i];   msg += ss.str();  ss.str("");   ss.clear();
                msg += "\n";
            }
            output->toplog(msg);
        }
        
        delete [] QkT;
        
        
        fMC_justreset = false; // this ensures that at the next time step the rate vectors in each subcomm are read and not compiled again
        
        if (me == MASTER) {
            QkTint += Qktot * Dt;
            fprintf(screen,"\n\n MASTER's cumulative rate integrated over time is: %f \n",QkTint);
        }
        
        
        
        
        
        
        // if a KMC process was chosen, master must tell all subcomms who the selected subcomm is and pass to them the cumrate value to check in the cumrate vector
        if (KMCexecute == 1){
            
            // all processors send their QkTC[end] to the QkTCall array in master with size nproc. All slaves should be sending a zero at this point, and submasters a value
            if (me>0){
                int dest = 0;
                MPI_Send(&QkTC[nKMCs-1], 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            }
            if (me == MASTER) {
                QkTCall[0] = QkTC[nKMCs-1];
                for (int source=1; source<universe->nprocs; source++) {
                    MPI_Recv(&QkTCall[source], 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD,&status);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            
            int chProc = -1;   // the ID of the processor (submaster) running the fix of the selected event
            double CRCred = 0.;  // same as CRC defined below but reduced to the processor containing the event to be executed
            int chScom = -1;    // color of Subcomm of chosen processor
            // master cumulates QkTCall
            if (me == MASTER) {
                QkTCallC[0] = QkTCall[0];
                fprintf(screen,"\n\n MASTER's cumulative rates from all procs: \npID  QkTCall   QkTCallC  \n%d   %f   %f\n",0,QkTCall[0],QkTCallC[0]);
                for (int i=1; i<universe->nprocs; i++)  {
                    QkTCallC[i] = QkTCallC[i-1] + QkTCall[i];
                    fprintf(screen,"%d   %f   %f\n",i,QkTCall[i],QkTCallC[i]);
                }
                
                // using a new random number the master picks an event in QkTCallC cumulated array and communicates to all procs the corresponding processor (a submaster) and the random number minus QkTallC[pos-1]
                double u2 = randm->pick01();
                double CRC = u2 * QkTCallC[universe->nprocs-1];  // cumulative rate integral of chosen event: used to identify event to execute
                QkTint = 0.;
                
                fprintf(screen,"\n MASTER chose rate %f \n",CRC);
                
                if (CRC < QkTCallC[0]) {
                    chProc = 0;
                    CRCred = CRC;
                }
                else {
                    bool chfound = false;
                    for (int i=1; i<universe->nprocs; i++)  {
                        if (!chfound && CRC < QkTCallC[i]) {
                            chProc = i;
                            chfound = true;
                            CRCred = CRC - QkTCallC[i-1];
                        }
                    }
                }
                chScom = universe->color_each[chProc];
                
                fprintf(screen,"\n MASTER sends reduced rate %f to proc %d, in subcomm %s \n",CRCred,chProc,universe->SCnames[chScom].c_str());

                for (int source=1; source<universe->nprocs; source++) {
                    MPI_Send(&chProc, 1, MPI_INT, source, 1, MPI_COMM_WORLD);
                    MPI_Send(&CRCred, 1, MPI_DOUBLE, source, 2, MPI_COMM_WORLD);
                    MPI_Send(&chScom, 1, MPI_INT, source, 3, MPI_COMM_WORLD);
                }
            }
            
            if (me>0){
                int dest = 0;
                MPI_Recv(&chProc, 1, MPI_INT, dest, 1, MPI_COMM_WORLD,&status);
                MPI_Recv(&CRCred, 1, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD,&status);
                MPI_Recv(&chScom, 1, MPI_INT, dest, 3, MPI_COMM_WORLD,&status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            
            int EVtype = -1;  // event type: will be 0 if delete, 1 if nucleate..
            int EVpID = -1; // ID of particle selected for deletion or nucleation
            int EVafixID = -1; // absolute fix ID of particle selected for deletion or nucleation
            double EVpDIAM = -1.; // Diameter of particle selected for nucleation
            double EVpXpos = -1.; // Position of particle selected for nucleation
            double EVpYpos = -1.;
            double EVpZpos = -1.;
            int EVpTYPE = -1;    // type of particle selected for nucleation
            
            // the chosen processor (should be a submaster) determines whether the chosen event is delete, nucleate, or else, and records id ot the particle to be deleted/created and global id of the fix calling the event (this latter contains the info to later update the solution)
            if (me==chProc) {
                if (CRCred < QkTC[0]) {
                    EVtype = 0;
                }
                else {
                    bool typefound = false;
                    for (int i=1; i<nKMCs; i++)  {
                        if (!typefound && CRCred < QkTC[i]) {
                            EVtype = i;
                            typefound = true;
                            CRCred = CRCred - QkTC[i-1];
                        }
                    }
                }
                
                fprintf(screen,"\n PROC %d has these cumulated rates QkTC: %f %f \n",me,QkTC[0],QkTC[1]);

                fprintf(screen,"\n PROC %d choses event type %d with reduced rate %f \n",me,EVtype,CRCred);

                
                if (EVtype==0) {
                    fix_del->findEV(CRCred);
                    EVpID = fix_del->EVpID;   //more general to do this than returning by reference
                    EVafixID = fix_del->EVafixID;
                    //fix_del->EVpID = -1;   //(cleaning it up for safety)
                    //fprintf(screen,"\n PROC %d cumulated and T-integrated rates are:\n",me);
                    //for (int i=0; i<fix_del->CratesT.size(); i++) {
                      //fprintf(screen,"%f\n",fix_del->CratesT[i]);
                    //}
                    //fprintf(screen,"\n PROC %d choses particle %d to delete using fix %s\n",me,EVpID,fix->afKMCname[EVafixID].c_str());
                }
                else if (EVtype==1) {
                    fix_nucl->findEV(CRCred);
                    EVpID = fix_nucl->EVpID;
                    EVafixID = fix_nucl->EVafixID;
                    //fix_nucl->EVpID = -1;   //(cleaning it up for safety)
                    //fprintf(screen,"\n PROC %d choses particle %d to nucleate using fix %s\n",me,EVpID,fix->afKMCname[EVafixID].c_str());
                }
            }
            
            
            // chosen processor (me == chProc) sends data to execute events to others, and non chosen processors (me != chProc) receive them
            if (me==chProc) {
                // chosen processor sends event type
                for (int dest=0; dest<universe->nprocs; dest++) {
                    if (dest!=me) {
                        MPI_Send(&EVtype, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
                    }
                }
                // depending on event type, chosen processor sends info
                if (EVtype==0 || EVtype==1) {   // deletion and nucleation event only pass particle and fix IDs for now
                    for (int dest=0; dest<universe->nprocs; dest++) {
                        if (dest!=me) {
                            MPI_Send(&EVpID, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
                            MPI_Send(&EVafixID, 1, MPI_INT, dest, 3, MPI_COMM_WORLD);
                        }
                    }
                }
            }
            else {
                // other processors receive evnt type
                int source = chProc;  //NB: master told everybody who the chProc was
                MPI_Recv(&EVtype, 1, MPI_INT, source, 1, MPI_COMM_WORLD,&status);
                if (EVtype==0 || EVtype==1) {
                        MPI_Recv(&EVpID, 1, MPI_INT, source, 2, MPI_COMM_WORLD,&status);
                        MPI_Recv(&EVafixID, 1, MPI_INT, source, 3, MPI_COMM_WORLD,&status);
                }
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            
            // Sucommunicator where event takes place extracts from lammps any other quantity and passes them to all processors for everyone to carry out the event (nothing else than the already known for deletion events, whereas nucleation requires position and shape parameters)
            if (universe->color == chScom) {
                if (EVtype==0) {// do nothing - deletion only needs particle ID and all processors already know EVpID
                }
                else if(EVtype==1){  // nucleation instead requires chosen particle position and shape, to be extracted from LAMMPS - all procs in subcomm call a fix_nucleate function doing this
                    // for now I am assuming spherical shape and thus associated parameters, but in the future might want to extend to ellipsoids at least..
                    fix_nucl->extractPROPS(EVpID,EVafixID);
                    //fprintf(screen,"\n PROC %d started with ptype = %d\n",me,EVpTYPE);
                    EVpDIAM = fix_nucl->EVpDIAM;
                    EVpXpos = fix_nucl->EVpXpos;
                    EVpYpos = fix_nucl->EVpYpos;
                    EVpZpos = fix_nucl->EVpZpos;
                    EVpTYPE = fix_nucl->EVpTYPE;
                    //fprintf(screen,"\n PROC %d has identified ptype = %d\n",me,EVpTYPE);
                }
            }
            
            // submaster with selected event passes additional particle info to all other processors (though this includes processors in same subcomm, which actually know the info already
            if(me ==chProc){
                for (int dest=0; dest<universe->nprocs; dest++) {
                    if (dest!=me) {
                        MPI_Send(&EVpDIAM, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
                        MPI_Send(&EVpXpos, 1, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
                        MPI_Send(&EVpYpos, 1, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
                        MPI_Send(&EVpZpos, 1, MPI_DOUBLE, dest, 4, MPI_COMM_WORLD);
                        MPI_Send(&EVpTYPE, 1, MPI_INT, dest, 5, MPI_COMM_WORLD);
                    }
                }
            }
            else{
                // other processors receive evnt type
                int source = chProc;  //NB: master told everybody who the chProc was
                MPI_Recv(&EVpDIAM, 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD,&status);
                MPI_Recv(&EVpXpos, 1, MPI_DOUBLE, source, 2, MPI_COMM_WORLD,&status);
                MPI_Recv(&EVpYpos, 1, MPI_DOUBLE, source, 3, MPI_COMM_WORLD,&status);
                MPI_Recv(&EVpZpos, 1, MPI_DOUBLE, source, 4, MPI_COMM_WORLD,&status);
                MPI_Recv(&EVpTYPE, 1, MPI_INT, source, 5, MPI_COMM_WORLD,&status);
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            /*
            if (EVtype==0) {
                fprintf(screen,"\n PROC %d has particle deletion parameters: EVpID = %d ; fix = %s\n",me,EVpID,fix->afKMCname[EVafixID].c_str());
            }
            else if  (EVtype==1){
                fprintf(screen,"\n PROC %d has particle nucleation parameters: EVpID = %d ; fix = %s ; diam = %f ; XYZ = %f %f %f ; ptype = %d\n",me,EVpID,fix->afKMCname[EVafixID].c_str(),EVpDIAM,EVpXpos,EVpYpos,EVpZpos,EVpTYPE);
            }
            */
            
            
            
            
            // executing the event - all processors
            if (EVtype==0){  // if a deletion..
                double pV;   //part vol: to be used to update the solution
                // processors with active lammps delete the selected particle
                if (lammpsIO->lammps_active) {
                    for (int i=0; i<fix->fKMCtype.size(); i++){
                        if (strcmp(fix->fKMCtype[i].c_str(),"nucleate")==0){
                            fix_nucl->delete_trial(i);  // delete trial particles in all fixes, so that all subcomms and processors will have same number of particles. Then run a lammps 0 just to make sure..
                        }
                    }
                    lammpsIO->lammpsdo("run 0");
                    fix_del->execute(EVpID);
                    pV = fix_del->partVol;
                    fprintf(screen,"\n PROC %d Volume of particle to be deleted: %f\n",me,pV);

                    // if current processor is submaster of first subcomm with lammps active, it sends around the particle volume
                    if (me==universe->flampID) {
                        for (int dest =0; dest<universe->nprocs; dest++) {
                            if (dest!=me) {
                                MPI_Send(&pV, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                
                if (me!=universe->flampID) {
                    int source = universe->flampID;
                    MPI_Recv(&pV, 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD,&status);
                }
                
                MPI_Barrier(MPI_COMM_WORLD);
                
                // all processors update the solution
                solution->update(EVafixID,pV,EVtype);
            }
            else if (EVtype==1){  // if a nucleations...
                double pV;   //part vol: to be used to update the solution
                // processors with active lammps nucleate the selected particle with assigned ID
                if (lammpsIO->lammps_active) {
                    for (int i=0; i<fix->fKMCtype.size(); i++){
                        if (strcmp(fix->fKMCtype[i].c_str(),"nucleate")==0){
                            fix_nucl->delete_trial(i);  // delete trial particles in all fixes, so that all subcomms and processors will have same number of particles. Then run a lammps 0 just to make sure..
                        }
                    }
                    lammpsIO->lammpsdo("run 0");
                    fix_nucl->execute(EVpID,EVafixID,EVpTYPE,EVpDIAM,EVpXpos,EVpYpos,EVpZpos);   // FOR SPHERICAL... another function with more arguments would be needed for ellipsoids
                    pV = fix_nucl->partVol;
                    fprintf(screen,"\n PROC %d Volume of particle to be nucleated: %e\n",me,pV);
                    
                    // if current processor is submaster of first subcomm with lammps active, it sends around the particle volume
                    if (me==universe->flampID) {
                        for (int dest =0; dest<universe->nprocs; dest++) {
                            if (dest!=me) {
                                MPI_Send(&pV, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                
                if (me!=universe->flampID) {
                    int source = universe->flampID;
                    MPI_Recv(&pV, 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD,&status);
                }
                
                MPI_Barrier(MPI_COMM_WORLD);
                
                // all processors update the solution
                fprintf(screen,"\n PROC %d (submaster) Updating the solution \n",me);
                solution->update(EVafixID,pV,EVtype);
                fprintf(screen,"\n PROC %d (submaster) solution updated\n",me);

            }
            
            // having carried out the discrete event, selcect random number to be used in determining next FOT of discrete events
            u = randm->pick01();
    
            // ********* NEED TO CREATE A RESET VEC FUCNTION SO THAT YOU CAN CALL IT HERE TOO *********
            // all processors zero their KMC-related rates and cumulated rates
            //if (universe->key ==0){}
            reset_event_vec = true;     //this is a harsh way of doing this. In future we should be more surgical and just update rate vectors etc. In this way everything is reset.
            // NB!!: if you do not reset the vectors after one KMC event is realised, you will have problems with IDs because the default option in delete_atoms is "compress yes", which reassigns IDs when a particle is deleted
        }
        delete [] QkTC;
        if (me==MASTER){
            delete [] QkTCall;
            delete [] QkTCallC;
        }
	
	if (Cexecute) {
	  // reset the position of trial atoms
	  for (int i=0; i<fix->fKMCtype.size(); i++) {
	    if (strcmp(fix->fKMCtype[i].c_str(),"nucleate")==0) {
	      fix_nucl->reset_trial_pos(i);
	    }
	  }
#ifdef MASKE_WITH_NUFEB
	  if (universe->color == SC2exec) {
	    if (type == "nufeb") {
	      fix_nufeb->execute(Cpid2exec, SC2exec);
	      fix->Cleval[Cpid2exec] = msk->tempo;
	    }
	  } else {
	    if (type == "nufeb") {
	      // Receive new concentrations from nufeb
	      for (int i = 0; i < chem->Nmol; i++) {
		if (chem->mol_nufeb[i] > 0) { // if points to a valid nufeb chemical species
		  double conc = 0;
		  MPI_Bcast(&conc, 1, MPI_DOUBLE, universe->subMS[SC2exec], MPI_COMM_WORLD);
		  chem->mol_cins[i] = conc;
		}
	      }
	    }
	  }
	  if (type == "nufeb") {
	    fix_nufeb->exchange(Cpid2exec, SC2exec);
	  }
#endif
	}

        if (endofrun==1){
            // delete all trial particles associated with nucleation fixes
            for (int i=0; i<fix->fKMCtype.size(); i++){
                if (strcmp(fix->fKMCtype[i].c_str(),"nucleate")==0){
                    fix_nucl->delete_trial(i);  // delete trial particles in all fixes, so that all subcomms and processors will have same number of particles. Then run a lammps 0 just to make sure..
                }
            }
            lammpsIO->lammpsdo("run 0");
       
            // execute all continuous processes up to now
            // To be implemeneted because to date there are no continuou processes implemented..
            // Not sure what to do with cumulated KMC rates till here: Loosing them is a pity/error, but keeping them is risky because the user can unfix and do things before the next Krun, leading to inconsistencies..
        }
        
        // number of just completed step
        msk->step = msk->step+1;
        
        #ifdef MASKE_WITH_SPECIATION
            // call speciation if at appropriate step
            for (int i = 0; i<spec->specID.size();i++){
                if (msk->step % spec->spec_every[i] == 0) {
                    spec->dospec(i);
                }
            }
        #endif
        
        
        // relax if at appropriate step
        for (int i = 0; i<relax->rlxID.size();i++){
            if (msk->step % relax->rlx_every[i] == 0) {
                relax->dorelax(i);
            }
        }

        // If a KMC event was not executed
        // (must stay here after relax, in case a continuous event was carried out and then the relax changed the box)
        if (KMCexecute == 0){
            // NUKK: Reset positions of trial particles in all fix nucleates to their initial lattice site, plus any distortion occurred to the simulation box since last accepted KMC event
            // A loop considering all fixes on processor and, if the fix is a nucleate, call a function within fix_nucleate to reset the poision of the trial particles duly
        }
        
        // write thermo output at appropriate step
        if (output->th_every>0 && msk->step % output->th_every == 0) {
	  output->writethermo();
        }
        // write dump output at appropriate step
        for (int i = 0; i<output->dumpID.size();i++){
            if (msk->step % output->dump_every[i] == 0) {
	        output->writedump(i);
                output->dump_first[i]=false;  // from not on, not the first time the dump is written anymore
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        //   ***** END OF LOOP

        }
    
    fprintf(screen,"\n Proc[%d]: KRUN loop completed succesfully. Moving on.. \n",me);

}









    
// ---------------------------------------------------------------
// Printing info about the class (possibly useful for debugging)
void Krun::printall()
{
	fprintf(screen,"\n---------ALL ABOUT KRUN----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
