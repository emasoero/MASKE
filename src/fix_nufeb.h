#ifndef FIX_NUFEB_H
#define FIX_NUFEB_H

#include "pointers.h"
//#include <string>
#include <vector>
#include "mpi.h"
//#include <unistd.h>   //just for the sleep() function
//#include <string>
//#include <vector>
/*#include <stdlib.h>
  #include <stdio.h>
  #include <iostream>
  #include <fstream>
*/

#define MASTER 0
#define M_PI 3.14159265358979323846

namespace MASKE_NS {
	
  class Fix_nufeb : protected Pointers {
  public:
    Fix_nufeb(class MASKE *);
    ~Fix_nufeb();
    
    void init(int pos);
    double getDT(int);       // function computing time increment for the continuum process
    void printall();
    void execute(int);

    int group; // bacteria group
  };
}

#endif
