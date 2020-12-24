# MASKE

The basic version of MASKE only requires the LAMMPS submodule.
If you want to include speciation in your simulation, than you need LAMMPS and PHREEQC.
If you want to have bacteria too, you need LAMMPS, PHREEQC, and NUFEB.

## Dependecies

| Name   | Description                        | Version                          |
|--------|------------------------------------|----------------------------------|
| Git    | Distributed version control system | Most recent versions should work |
| CMake  | Cross-platform build tool          | 3.0 or above                     |
| MPI    | Implementation of the MPI standard | Same requirements as LAMMPS      |

## Download

### Downloading basic MASKE (only LAMMPS, no PHREEQC nor NUFEB)

Use:
```
$ cd <path> 
$ git clone https://github.com/emasoero/MASKE
$ cd MASKE
$ git submodule update --init lammps
```
*path* should be the place where you want the source code to be downloaded.

### Downloading MASKE (and LAMMPS) with PHREEQC (no NUFEB)

Use:

NO IDEA HERE.... DENIS?

### Downloading MASKE with NUFEB (and LAMMPS and PHREEQC)

Use:
```
$ cd <path> 
$ git clone --recurse-submodules https://github.com/emasoero/MASKE
```

## Folder layout

After download is completed the MASKE folder should contain the following folders:

```
$ cd MASKE 
$ ls
```

| Folder | Description                |
|--------|----------------------------|
| src    | Source files               |
| lib    | Library dependancies files |
| tests  | Collection of test cases   |
| bin    | Built executables          |
| lammps | LAMMPS source folder       |

DENIS: there is no lib, nor bin, but there is NUFEB.... TO BE CORRECTED?

## Building

We rely on [CMake](https://cmake.org) to provide cross-platform build support.

### Linux and MacOS

First we compile LAMMPS. From within the MASKE folder, use:
```
$ cd lammps
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DWITH_JPEG=Off -DWITH_PNG=Off -DWITH_FFMPEG=Off -DBUILD_LIB=On -DBUILD_OMP=Off -DPKG_MISC=yes -DPKG_USER-MASKE=yes ../cmake
$ make -j
$ cd ../..
$ mkdir build
```

If you want to install NUFEB, disregard the box below and move to the section "Enabling NUFEB" below.

If you want to install MASKE without NUFEB, Use:
```
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make -j
```
If you want to enable speciation by PHREEQC, add the cmake option  -DCMAKE_WITH_SPECIATION=On

DENIS: is this correct about PHREEQC? What should the user do to actually make sure they have a working PHREEQC?

### Newcastle University Rocket cluster

You will need first to load the required modules:
```
$ module load CMake
$ module load OpenMPI
```
Then follow the Linux build instructions.

### Enabling NUFEB

MASKE can be coupled with microbial simulation provided by LAMMPS's [NUFEB](https://github.com/nufeb/NUFEB) module. Optionally, it is possible to enable NUFEB with VTK support to output .vti and .vtu files for grid based data and per-atom data. You can use those files to visualise the simulation using [Paraview](https://www.paraview.org/). For more information on how to generate those files plase refer to the NUFEB user manual.

Use the following build steps:
```
$ cd NUFEB
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DWITH_NUFEB=On ..
```
You will get an CMake error after this but just ignore it.
```
$ cd ../NUFEB
```
If you don't want to enable VTK support do:
```
$ ./install.sh --enable-misc --static
$ cd ../build
$ cmake -DCMAKE_BUILD_TYPE=Release -DWITH_NUFEB=On ..
```
Otherwise do:
```
$ ./install.sh --enable-misc --enable-vtk --static
$ cd ../build
$ cmake -DCMAKE_BUILD_TYPE=Release -DWITH_NUFEB=On -DWITH_VTK=On -DVTK_DIR=$PWD/../NUFEB/thirdparty/vtk/vtk-build ..
```
Finally, compile MASKE using:
```
$ make -j
```

### Troubleshooting

- CMake can't find MPI compiler.

In the build folder type the following command:
```
$ cmake -DMPI_CXX_COMPILER=<path> ..
$ make -j
```
with *path* being the path to the desired MPI compiler.

- Compiler complains about c++11 standard.

If you are trying to compile LAMMPS when you get the error, in the lammps/build folder type the following command:
```
$ cmake -DCMAKE_CXX_FLAGS="-std=c++11" ../cmake
$ make -j
```

If you are trying to compile MASKE, in the MASKE/build folder type the following command:
```
$ cmake -DCMAKE_CXX_FLAGS="-std=c++11" ..
$ make -j
```

- On MacOS, ``make -j`` can sometimes create problems. Try using only ``make``, without the ``-j`` option

## Running

Let's try running an example located at the tests folder:
```
$ cd tests/test12-Lay_Nucl_2Grp_Sm_Mono_1pr
$ mpiexec -np 4 ../../build/maske input_read.dat 
```

### Running on Newcastle University Rocket cluster

Follow the instructions on the [Rocket webpage](https://services.ncl.ac.uk/itservice/research/hpc/). For a quick interactive run do:
```
$ cd tests/test12-Lay_Nucl_2Grp_Sm_Mono_1pr
$ srun -p interactive -n 4 --pty /bin/bash
$ mpiexec -np 4 ../../build/maske input_read.dat 
```
