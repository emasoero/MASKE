DIR = ./
PROGRAMS = $(DIR)main
CC = /opt/local/bin/mpicxx-mpich-mp
DIR       = OBJECTS/
DIREX     = RUN/


CFLAGS=-O3
LINKER= $(CC) $(CFLAGS)
 
LIBS=    -L/usr/lib -L/opt/local/lib -L/usr/local/lib  -L/opt/local/lib/mpich-mp -L/Users/emasoero/Documents/Programs/lammps-22Aug18-maske/src    -L/Users/nem42/Documents/Research/PROGRAMS/lammps-22Aug18-MASKE/src -L/Users/kumarancoopamootoo/Documents/lammps-11Aug17-maske/src -llammps_mac_mpi  -lmpich -lmpl -lpthread -lstdc++ 


#additional lammps libraries that you might add at some point if you want to use extra packages (need recompiling lammps as library though)
# -lfftw

INCLUDES=  -I/usr/include -I/opt/local/include -I/usr/local/include -I/opt/local/include/mpich-mp -I/Users/emasoero/Documents/Programs/lammps-22Aug18-maske/src  -I/Users/nem42/Documents/Research/PROGRAMS/lammps-22Aug18-MASKE/src -I/Users/kumarancoopamootoo/Documents/lammps-11aug17-maske/src


# all: $(PROGRAMS)
# .cpp: ;  $(CC) $(CFLAGS) $(INCLUDES) $@.cpp $(LIBS) -o ../m3d

OBJ = $(DIR)main.o\
$(DIR)maske.o\
$(DIR)memory.o\
$(DIR)error.o\
$(DIR)inputmsk.o\
$(DIR)interactions.o\
$(DIR)universe.o\
$(DIR)DTnucleate.o\
$(DIR)lammpsIO.o\
$(DIR)chemistry.o\
$(DIR)particles.o\
$(DIR)simbox.o\
$(DIR)solution.o\
$(DIR)fix.o\
$(DIR)krun.o\
$(DIR)randm.o\
$(DIR)fix_delete.o\
$(DIR)output.o\
$(DIR)fix_Cfoo.o\
$(DIR)relax.o\
$(DIR)fix_nucleate.o\
$(DIR)block.o\
$(DIR)store.o

MAKEFILE = makefile

$(DIR)%.o: %.cpp $(MAKEFILE)
	$(CC) -c $(CFLAGS) -o $@ $(INCLUDES) $*.cpp

$(DIREX)maske:  $(OBJ) $(MAKEFILE)
	$(LINKER) $(OBJ)  -o $@ $(LIBS)

clean:  
	rm -f $(DIR)*.o	
