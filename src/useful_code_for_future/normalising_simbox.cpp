#include "simbox.h"
//#include "universe.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "clinker.h"
#include "pyramids.h"
#include "hydration.h"
#include "memory.h"
#include "csh.h"
#include "c2s.h"
#include "c3s.h"
#include "water.h"
#include <string.h>
*/

using namespace MASKE_NS;


Simbox::Simbox(MASKE *maske) : Pointers(maske)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    /*MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating lammpsIO class\n");

    lammps_active=false;
    units = "real";
    atomstyle = "ellipsoid";
     */
    for(int i=0;i<3;i++) boundary.push_back("p");
   
    xyzlo[0]=0.;
    xyzlo[1]=0.;
    xyzlo[2]=0.;
    
    xyzhi[0]=1.;
    xyzhi[1]=1.;
    xyzhi[2]=1.;
    
    xyztri[0]=0.;
    xyztri[1]=0.;
    xyztri[2]=0.;
    
    triclinic = false;
    
    
    mtri = (double**) malloc(3 * sizeof(double*));
    for (int i = 0; i < 3; i++) mtri[i] = (double*) malloc(3 * sizeof(double));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mtri[i][j] = 0;
            if (i == j) mtri[i][j] = 1.;
        }
    }
    
    mtriINV = (double**) malloc(3 * sizeof(double*));
    for (int i = 0; i < 3; i++) mtriINV[i] = (double*) malloc(3 * sizeof(double));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mtriINV[i][j] = 0;
            if (i == j) mtriINV[i][j] = 1.;
        }
    }
    

}



// ---------------------------------------------------------------
// Class destructor
Simbox::~Simbox()
{
    //delete lmp;
}





// ---------------------------------------------------------------
// Class destructor
void Simbox::computeTM()
{
    // calculate the transformation matrix from unit box to triclinic one, and the inverse from triclinic to unit
    mtri[0][0] = xyzhi[0] - xyzlo[0];
    mtri[0][1] = xyztri[0];
    mtri[1][1] = xyzhi[1] - xyzlo[1];
    mtri[0][2] = xyztri[1];
    mtri[1][2] = xyztri[2];
    mtri[2][2] = xyzhi[2] - xyzlo[2];
    
    //fprintf(screen, "\n\n The normalized --> triclinic transformation matrix is \n %f %f %f \n %f %f %f \n %f %f %f \n\n", mtri[0][0], mtri[0][1], mtri[0][2], mtri[1][0], mtri[1][1], mtri[1][2], mtri[2][0], mtri[2][1], mtri[2][2]);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mtriINV[i][j] = mtri[i][j];
        }
    }
    
    MatInv(mtriINV);
    
    //fprintf(screen, "\n\n The triclinic --> normalized transformation matrix is \n %f %f %f \n %f %f %f \n %f %f %f \n\n", mtriINV[0][0], mtriINV[0][1], mtriINV[0][2], mtriINV[1][0], mtriINV[1][1], mtriINV[1][2], mtriINV[2][0], mtriINV[2][1], mtriINV[2][2]);
    
    
    // test the triclinic to 111 normalized and inverse vector transformations
    double Vtest[3] = {1, 1, 1};
    N2Ttrans(mtri, Vtest, xyzlo);
    //fprintf(screen, "\n\n The 111 nomralized vector, in the triclinic box becomes \n %f %f %f \n\n", Vtest[0], Vtest[1], Vtest[2]);
    T2Ntrans(mtriINV, Vtest, xyzlo);
    //fprintf(screen, "\n\n And now the 111 vector should be recovered \n %f %f %f \n\n", Vtest[0], Vtest[1], Vtest[2]);

}





// ---------------------------------------------------------------
// function to make a matrix x column vector product
void Simbox::MVprod(double **M, double *v)
{
    double vt[3];
    for (int i = 0; i < 3; i++) {
        vt[i] = M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2];
    }
    for (int i = 0; i < 3; i++) {
        v[i] = vt[i];
    }
    return;
}


// ---------------------------------------------------------------
// function to transform a vector from a triclinic box with origin in xyzlo to a vector in normalized 111 orthogonal box
void Simbox::N2Ttrans(double **M, double *v, double *ori)
{
    MVprod(M, v);
    for (int i = 0; i < 3; i++) {
        v[i] = v[i] + ori[i];
    }
    return;
}

// ---------------------------------------------------------------
// function to transform a vector from a normalized 111 orthogonal to a triclinic box with origin in xyzlo
void Simbox::T2Ntrans(double **M, double *v, double *ori)
{
    for (int i = 0; i < 3; i++) {
        v[i] = v[i] - ori[i];
    }
    MVprod(M, v);
    return;
}


// ---------------------------------------------------------------
void Simbox::MatInv(double **M)
{
    //PAX=LUX=I algorithm
    double a[3][3];
    double d[3], sw[3];
    double temp[3][3], P[3][3];
    double ident[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double zeros[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            temp[i][j] = M[i][j];
        }
    }
    
    for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
            P[ii][jj] = ident[ii][jj];
        }
    }
    
    
    //swap first line if it starts with zero. Record correspondent permutation in P
    if (temp[0][0] == 0.) {
        if (temp[1][0] != 0.) {
            sw[0] = temp[1][0]; sw[1] = temp[1][1]; sw[2] = temp[1][2];
            temp[1][0] = temp[0][0]; temp[1][1] = temp[0][1]; temp[1][2] = temp[0][2];
            temp[0][0] = sw[0]; temp[0][1] = sw[1]; temp[0][2] = sw[2];
            P[0][0] = 0.; P[0][1] = 1.; P[1][0] = 1.; P[1][1] = 0.;
        }
        else {
            sw[0] = temp[2][0]; sw[1] = temp[2][1]; sw[2] = temp[2][2];
            temp[2][0] = temp[0][0]; temp[2][1] = temp[0][1]; temp[2][2] = temp[0][2];
            temp[0][0] = sw[0]; temp[0][1] = sw[1]; temp[0][2] = sw[2];
            P[0][0] = 0.; P[0][2] = 1.; P[2][0] = 1.; P[2][2] = 0.;
        }
    }
    
    //swap second and third line if necessary, and update the permutation
    if (temp[1][1] - temp[1][0] / temp[0][0]*temp[0][1] == 0.) {
        double P2[3][3];
        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                P2[ii][jj] = zeros[ii][jj];
            }
        }
        sw[0] = temp[2][0]; sw[1] = temp[2][1]; sw[2] = temp[2][2];
        temp[2][0] = temp[1][0]; temp[2][1] = temp[1][1]; temp[2][2] = temp[1][2];
        temp[1][0] = sw[0]; temp[1][1] = sw[1]; temp[1][2] = sw[2];
        P2[0][0] = 1.; P2[1][2] = 1.; P2[2][1] = 1.;
        
        double Ptemp[3][3];
        for (int lp = 0; lp < 3; lp++) {
            for (int k = 0; k < 3; k++) {
                Ptemp[lp][k] = P2[lp][0] * P[0][k] + P2[lp][1] * P[1][k] + P2[lp][2] * P[2][k];
            }
        }
        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                P[ii][jj] = Ptemp[ii][jj];
            }
        }
    }
    
    temp[1][0] /= temp[0][0];
    temp[1][1] -= temp[0][1] * temp[1][0];
    temp[1][2] -= temp[0][2] * temp[1][0];
    
    temp[2][0] /= temp[0][0];
    temp[2][1] -= temp[0][1] * temp[2][0];
    temp[2][2] -= temp[0][2] * temp[2][0];
    
    temp[2][1] /= temp[1][1];
    temp[2][2] -= temp[1][2] * temp[2][1];
    
    for (int s = 0; s < 3; s++) {
        d[0] = ident[s][0];
        d[1] = ident[s][1] - temp[1][0] * d[0];
        d[2] = ident[s][2] - temp[2][0] * d[0] - temp[2][1] * d[1];
        
        a[2][s] = d[2] / temp[2][2];
        a[1][s] = (d[1] - temp[1][2] * a[2][s]) / temp[1][1];
        a[0][s] = (d[0] - temp[0][1] * a[1][s] - temp[0][2] * a[2][s]) / temp[0][0];
    }
    
    double atemp[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int lp = 0; lp < 3; lp++) {
        for (int k = 0; k < 3; k++) {
            atemp[lp][k] = a[lp][0] * P[0][k] + a[lp][1] * P[1][k] + a[lp][2] * P[2][k];
        }
    }
    for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
            M[ii][jj] = atemp[ii][jj];
        }
    }
    
    
    return;
}






// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Simbox::printall()
{
	fprintf(screen,"\n---------ALL ABOUT SIMBOX----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
