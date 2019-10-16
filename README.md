# MASKE

## Dependecies

| Name   | Description                        | Version                          |
|--------|------------------------------------|----------------------------------|
| Git    | Distributed version control system | Most recent versions should work |
| CMake  | Cross-platform build tool          | 3.0 or above                     |
| MPI    | Implementation of the MPI standard | Same requirements as LAMMPS      |
| JPEG   | Image library                      | Most recent versions should work |
| PNG    | Image library                      | Most recent versions should work |

## Download

Download MASKE source code using Git:
```
$ cd <path> 
$ git clone --recursive-submodules https://github.com/emasoero/MASKE
```
<path> should be the place where you want the source code downloaded.

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

## Building

We rely on [CMake](https://cmake.org) to provide cross-platform build support.

### Linux

```
$ cd lammps
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_LIB=on -DBUILD_OMP=off ../cmake
$ make -j
$ cd ../../
$ mkdir build
$ cd build
$ cmake ..
$ make -j
```
