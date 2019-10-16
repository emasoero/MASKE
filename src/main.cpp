#include "maske.h"
#include "mpi.h"
#include "inputmsk.h"
#include "universe.h"
#include <unistd.h>   //just for the sleep() function
//#include "clinker.h"
//#include "pyramids.h"
//
//#include "hydration.h"
//#include "air.h"
//#include "csh.h"


using namespace MASKE_NS;

int main(int argc, char **argv)
{
    
    MPI_Init(&argc,&argv);

    
	// Create an object from the fundamental (this object is essentially your whole simulation)
    MASKE *maske = new MASKE(argc,argv);
    
    
    
    // Read the input file (or, automatically, the restart file if there is one that is active)
    maske->inputmsk->file();
    

    
    // Allocate all the transition and continuum processes to the appropriate subcommunicator
    //maske->universe->populate();
    
    
    
    
    // Run the main loop
    //maske->mainloop();
    
    
    
    
    
    
    /*
     
	//read the inputcprs file and record values and flags
	setting->inputcprs->file();
	
	//import the clinker file
	setting->clinker->file();
	
	
	for (int i=0; i<setting->clinker->nclink_bins; i++) {
		fprintf(setting->screen,"clinker %f %f \n",setting->clinker->Ninbin[i],setting->clinker->diam[0][i]);
	}
	
	//compute initial surface and volume of clinker
	setting->clinker->init_surf_vol();
	
	//compute volume of the different phases in clinker, and check that they cover 100% of it
    //assign the specific volume of the clinker based on its composition
	setting->clinker->init_vol_phases();

	//generate a capillary pore distribution
	setting->pyramids->gen_pyramids();

	fprintf(setting->screen,"\n\nCOMPUTING INITIAL SURFACE AND WC\n ");
	//compute initial surface and volume of the capillary porosity, and so the initial volume of water
	setting->pyramids->init_surf_vol_wc();

	fprintf(setting->screen,"\n\nMESHING PYRAMIDS \n ");
	//mesh the pores with concentric spheres
	setting->pyramids->mesh();

	
	fprintf(setting->screen,"SETTING RZ\n\n ");
	fflush(setting->screen);
	//setting the random distribution of RZ thickness
	setting->pyramids->setRZ();
	
	fprintf(setting->screen,"SETTING NUCLETION RATES\n\n ");
	fflush(setting->screen);
	//setting the random distribution of RZ thickness
	setting->csh->setIb();

	//for (int i=0; i<setting->clinker->nclink_bins; i++) {
	//	fprintf(setting->screen,"clinker %i %f \n",setting->clinker->Ninbin[i],setting->clinker->diam[1][i]);
	//}
	
	
	fprintf(setting->screen,"\n\nSTARTING HYDRATION\n ");
	//initialize hydration
	setting->hydration->init();
	
	
	
	
	
	//for (int i=0; i<setting->clinker->nclink_bins; i++) {
	//	fprintf(setting->screen,"clinker %i %f \n",setting->clinker->Ninbin[i],setting->clinker->diam[2][i]);
	//}
	
	
	
	//advance hydration steps and compute everything
	int conta=0;
	fprintf(setting->screen,"Hydration completed at:\n ");
	for (int i=1; i< ((setting->hydration->Ndiss)+1); i++) {
		if ((double)i/(double)(setting->hydration->Ndiss) > (double)conta * 0.1) {
			conta++;
			fprintf(setting->screen,"%i%% \t", (int)((double)i/(double)(setting->hydration->Ndiss)*100.));
			fflush(setting->screen);
			fprintf(setting->screen,"%e \t %e \n", setting->air->vol,setting->clinker->vol);
		}
		setting->hydration->step(i);
	}
	
	*/
    
    
 

    
    
    
	maske -> printall();




    
    
    // Switch off the restart file
    //cipresso->restart->restart_off()
    
    //MPI_Comm_free(&(maske->universe->subcomm));
    
   
    
    delete maske;
    
    /*int sleeptime = 1;
    sleep(sleeptime);
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(maske->screen,"\n UE  PROC %d: eddig jo",maske->universe->me);
    //fprintf(screen,"\n PROC %d: eddig jo",me);
    sleep(sleeptime);
    MPI_Barrier(MPI_COMM_WORLD);
*/
        
    
    //MPI_Barrier(MPI_COMM_WORLD);


    
    //MPI_Comm_free(&(maske->universe->subcomm));
    MPI_Finalize();
   
    
     
	//return 0;
}

