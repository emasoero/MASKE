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
  src/block.cpp
  src/chemistry.cpp
  src/DTnucleate.cpp
  src/error.cpp
  src/fix_Cfoo.cpp
  src/fix.cpp
  src/fix_delete.cpp
  src/fix_nucleate.cpp
  src/inputmsk.cpp
  src/interactions.cpp
  src/krun.cpp
  src/lammpsIO.cpp
  src/main.cpp
  src/maske.cpp
  src/memory.cpp
  src/output.cpp
  src/particles.cpp
  src/randm.cpp
  src/relax.cpp
  src/simbox.cpp
  src/solution.cpp
  src/store.cpp
  src/universe.cpp)

if (WITH_NUFEB)
  file(COPY
    ${CMAKE_SOURCE_DIR}/lammps/src/min.cpp
    ${CMAKE_SOURCE_DIR}/lammps/src/min.h
    ${CMAKE_SOURCE_DIR}/lammps/src/USER-MASKE/min_maske.h
    ${CMAKE_SOURCE_DIR}/lammps/src/USER-MASKE/min_maske.cpp
    DESTINATION ${CMAKE_SOURCE_DIR}/NUFEB/lammps/src)
  find_library(LAMMPS_LIBRARIES lammps PATHS ${CMAKE_SOURCE_DIR}/NUFEB/lammps/src NO_DEFAULT_PATH)
  list(APPEND SOURCES src/fix_nufeb.cpp)
else()
  find_library(LAMMPS_LIBRARIES lammps PATHS ${CMAKE_SOURCE_DIR}/lammps/build NO_DEFAULT_PATH)
endif()

add_executable(maske ${SOURCES})
target_include_directories(maske PUBLIC ${MPI_CXX_INCLUDE_PATH} lammps/src)
target_link_libraries(maske ${MPI_CXX_LIBRARIES} ${LAMMPS_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${VTK_LIBRARIES})

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
