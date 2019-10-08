#ifndef MEMORY_H
#define MEMORY_H

#include "pointers.h"
//#include "stdio.h"
#include "stdlib.h"

namespace MASKE_NS {

class Memory : protected Pointers {
 public:
	Memory(class MASKE *);

	
	
	
	/* ----------------------------------------------------------------------
	 create a 2d array 
	 ------------------------------------------------------------------------- */
	
	template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2) 
    {
		int nbytes = ((int) sizeof(TYPE)) * n1*n2;
		TYPE *data = (TYPE *) malloc(nbytes);
		nbytes = ((int) sizeof(TYPE *)) * n1;
		array = (TYPE **) malloc(nbytes);
		
		int n = 0;
		for (int i = 0; i < n1; i++) {
			array[i] = &data[n];
			n += n2;
		}
		return array;
    }


	
	
	
	/* ----------------------------------------------------------------------
	 create a 3d array 
	 ------------------------------------------------------------------------- */
	
	template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, int n3) 
    {
		int nbytes = ((int) sizeof(TYPE)) * n1*n2*n3;
		TYPE *data = (TYPE *) malloc(nbytes);
		nbytes = ((int) sizeof(TYPE *)) * n1*n2;
		TYPE **plane = (TYPE **) malloc(nbytes);
		nbytes = ((int) sizeof(TYPE **)) * n1;
		array = (TYPE ***) malloc(nbytes);
		
		int i,j;
		int m;
		int n = 0;
		for (i = 0; i < n1; i++) {
			m = ((int) i) * n2;
			array[i] = &plane[m];
			for (j = 0; j < n2; j++) {
				plane[m+j] = &data[n];
				n += n3;
			}
		}
		return array;
    }

	
	
	
	/* ----------------------------------------------------------------------
	 create a 4d array 
	 ------------------------------------------------------------------------- */
	
	template <typename TYPE>
    TYPE ****create(TYPE ****&array, int n1, int n2, int n3, int n4)
    {
		int nbytes = ((int) sizeof(TYPE)) * n1*n2*n3*n4;
		TYPE *data = (TYPE *) malloc(nbytes);
		nbytes = ((int) sizeof(TYPE *)) * n1*n2*n3;
		TYPE **cube = (TYPE **) malloc(nbytes);
		nbytes = ((int) sizeof(TYPE **)) * n1*n2;
		TYPE ***plane = (TYPE ***) malloc(nbytes);
		nbytes = ((int) sizeof(TYPE ***)) * n1;
		array = (TYPE ****) malloc(nbytes);
		
		int i,j,k;
		int m1,m2,m3;
		int n = 0;
		for (i = 0; i < n1; i++) {
			m2 = ((int) i) * n2;
			array[i] = &plane[m2];
			for (j = 0; j < n2; j++) {
				m1 = ((int) i) * n2 + j;
				m2 = ((int) i) * n2*n3 + j*n3;
				plane[m1] = &cube[m2];
				for (k = 0; k < n3; k++) {
					m1 = ((int) i) * n2*n3 + j*n3 + k;
					cube[m1] = &data[n];
					n += n4;
				}
			}
		}
		return array;
    }
	
	
};
}

#endif
