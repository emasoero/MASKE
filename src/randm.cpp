#include "randm.h"
#include "universe.h"

using namespace MASKE_NS;

// ---------------------------------------------------------------
// Initialize class
Randm::Randm(MASKE *maske) : Pointers(maske)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    
    if (me == MASTER) fprintf(screen,"Generating randm class\n");
}




// ---------------------------------------------------------------
// Class destructor
Randm::~Randm()
{
    
}



// ---------------------------------------------------------------
// Seeds the random number generator for all subcomms
void Randm::seedit(int seme)
{
    
    srand(seme);
    MPI_Barrier(MPI_COMM_WORLD);

    /*int sleeptime;
    sleeptime = 2*me;
    sleep(sleeptime);
    fprintf(screen,"\n\n Processor %d, part of subcomm %d, has seme %d and first random number is %d \n",me,universe->color,seme,rand());
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
*/
    
    /*
    // master seeds its random generator and starts to communicate random numbers to other subcomm masters
    if (me==MASTER) {
        srand(seme);
        
        // array of random seeds to communicate
        double seed2scom [universe->nsc];
        // a random seed for each subcomm, including the MASTER's one
        for (int i=0; i<universe->nsc; i++) {
            seed2scom[i] = rand();
        }
        
        // maseter sends seed to each subcomm (all processor in same subcomm, i.e. same color, have same seed)
        for (int dest=1; dest<nprocs; dest++) {
            seme = seed2scom[universe->color_each[dest]];
            fprintf(screen,"\n\n MASTER sending seme %d to processor %d, part of subcom %d \n",seme,dest,universe->color_each[dest]);

            MPI_Send(&seme, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
        }
        
        seme = seed2scom[0];
    }
    
    if (me > MASTER) {
        MPI_Recv(&seme, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // all processors initialise their random generator. All processors in same subcomm have same seed, thus initialise their random generator equally
    srand(seme);
    
    int sleeptime;
    sleeptime = 2*me;
    sleep(sleeptime);
    fprintf(screen,"\n\n Processor %d, part of subcomm %d, has first random number is %d \n",me,universe->color,rand());
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
     */
    
}


// ---------------------------------------------------------------
// Picks a random number
double Randm::pick01()
{
    double r = ((double) rand() / (RAND_MAX));
    return r;
}




// ---------------------------------------------------------------
// Printing info (possibly useful for debugging)
void Randm::printall()
{
    fprintf(screen,"\n---------ALL ABOUT RANDM----------\n");
    //fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
    fprintf(screen,"---------------------------------------\n\n");
    
}
