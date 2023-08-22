cmake_minimum_required(VERSION 3.0)

project(MASKE VERSION 1.0 LANGUAGES CXX)

find_package(MPI REQUIRED)
if (MPI_FOUND)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()
if (WITH_VTK)
  find_package(VTK)
endif()

set(SOURCES 
  src/chemistry.cpp
  src/error.cpp
  src/fix.cpp
  src/fix_CmstoreLMP.cpp
  src/fix_EmstoreLMP.cpp
  src/fix_delete.cpp
  src/fix_nucleate.cpp
  src/inputmsk.cpp
  src/krun.cpp
  src/lammpsIO.cpp
  src/main.cpp
  src/maske.cpp
  src/memory.cpp
  src/output.cpp
  src/randm.cpp
  src/relax.cpp
  src/setconc.cpp
  src/solution.cpp
  src/store.cpp
  src/universe.cpp)

if (WITH_NUFEB)
  file(COPY
    ${CMAKE_SOURCE_DIR}/lammps/src/min.cpp
    ${CMAKE_SOURCE_DIR}/lammps/src/min.h
    ${CMAKE_SOURCE_DIR}/lammps/src/USER-MASKE/min_maske.h
    ${CMAKE_SOURCE_DIR}/lammps/src/USER-MASKE/min_maske.cpp
    ${CMAKE_SOURCE_DIR}/lammps/src/USER-MASKE/min_quickmaske.h
    ${CMAKE_SOURCE_DIR}/lammps/src/USER-MASKE/min_quickmaske.cpp
    DESTINATION ${CMAKE_SOURCE_DIR}/NUFEB/lammps/src)
  find_library(LAMMPS_LIBRARIES lammps PATHS ${CMAKE_SOURCE_DIR}/NUFEB/lammps/src NO_DEFAULT_PATH)
  list(APPEND SOURCES src/fix_nufeb.cpp)
else()
  find_library(LAMMPS_LIBRARIES lammps PATHS ${CMAKE_SOURCE_DIR}/lammps/build NO_DEFAULT_PATH)
endif()

if (WITH_SPECIATION)
  find_library(PHREEQC_LIBRARIES iphreeqc PATHS ${CMAKE_SOURCE_DIR}/lib/iphreeqc/lib NO_DEFAULT_PATH)
  list(APPEND SOURCES src/spec.cpp)
endif()

add_executable(maske ${SOURCES})
target_include_directories(maske PUBLIC ${MPI_CXX_INCLUDE_PATH} lammps/src)
target_link_libraries(maske ${MPI_CXX_LIBRARIES} ${LAMMPS_LIBRARIES} ${PHREEQC_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${VTK_LIBRARIES})

if (MASKE_WAIT_ATTACH)
  target_compile_definitions(maske PUBLIC MASKE_WAIT_ATTACH)
endif()

if (MASKE_WAIT_ATTACH_RANK GREATER_EQUAL 0)
  target_compile_definitions(maske PUBLIC MASKE_WAIT_ATTACH_RANK=${MASKE_WAIT_ATTACH_RANK})
endif()

if (WITH_NUFEB)
  target_include_directories(maske PUBLIC NUFEB/src/USER-NUFEB)
  target_compile_definitions(maske PUBLIC MASKE_WITH_NUFEB)
endif()

if (WITH_SPECIATION)
  target_include_directories(maske PUBLIC ${CMAKE_SOURCE_DIR}/lib/iphreeqc-3.6.2-15100/include)
  target_compile_definitions(maske PUBLIC MASKE_WITH_SPECIATION)
endif()
