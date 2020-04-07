#ifndef SET_POINTERS_H
#define SET_POINTERS_H

#include "maske.h"

namespace MASKE_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define PI 3.1415926535897932384626433832795

class Pointers {

    
public:

    Pointers(MASKE *ptr) :
    msk(ptr),
    memory(ptr->memory),
    error(ptr->error),
    inputmsk(ptr->inputmsk),
    interact(ptr->interact),
    universe(ptr->universe),
    nucleate(ptr->nucleate),
    lammpsIO(ptr->lammpsIO),
    chem(ptr->chem),
    particles(ptr->particles),
    simbox(ptr->simbox),
    solution(ptr->solution),
    fix(ptr->fix),
    fix_del(ptr->fix_del),
#ifdef MASKE_WITH_NUFEB
    fix_nufeb(ptr->fix_nufeb),
#endif
    krun(ptr->krun),
    randm(ptr->randm),
    output(ptr->output),
    screen(ptr->screen),
    plog(ptr->plog),
    thermo(ptr->thermo),
    fix_cfoo(ptr->fix_cfoo),
    relax(ptr->relax),
    block(ptr->block),
    store(ptr->store),
    fix_nucl(ptr->fix_nucl) {}
        
    virtual ~Pointers() {}
        
protected:
    MASKE *msk;
    Memory *&memory;
    Error *&error;
    Inputmsk *&inputmsk;
    Universe *&universe;
    DTnucleate *&nucleate;
    LammpsIO *&lammpsIO;
    Chemistry *&chem;
    Particles *&particles;
    Simbox *&simbox;
    Solution *&solution;
    Interactions *&interact;
    Fix *&fix;
    Fix_delete *&fix_del;
    Fix_Cfoo *&fix_cfoo;
#ifdef MASKE_WITH_NUFEB
    Fix_nufeb *&fix_nufeb;
#endif
    Relax *&relax;
    Krun *&krun;
    Randm *&randm;
    Output *&output;
    Fix_nucleate *&fix_nucl;
    Block *&block;
    Store *&store;
    FILE *&screen;
    FILE *&plog;
    FILE *&thermo;
};
     
}

#endif
