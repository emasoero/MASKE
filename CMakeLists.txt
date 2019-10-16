cmake_minimum_required(VERSION 3.0)

project(MASKE VERSION 1.0 LANGUAGES CXX)

find_package(MPI REQUIRED)
if (MPI_FOUND)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()
find_package(JPEG REQUIRED)
find_package(PNG REQUIRED)

find_library(LAMMPS_LIBRARIES lammps ${CMAKE_SOURCE_DIR}/lammps/build)

add_executable(maske
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
target_include_directories(maske PUBLIC ${MPI_CXX_INCLUDE_PATH} lammps/src)
target_link_libraries(maske ${MPI_CXX_LIBRARIES} ${LAMMPS_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES})
