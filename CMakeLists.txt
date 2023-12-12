cmake_minimum_required(VERSION 3.18...3.28)
# Must enable languages here so that calls after FetchContent works
project(AplusC LANGUAGES C CXX CUDA)

include(CTest)
include(FetchContent)

# -- Externals
# Set up default options here (TBD if we want users to set these)
# Or.. in CMake setting files
set(CELERITAS_USE_CUDA ON)
set(CELERITAS_USE_Geant4 ON)
set(CELERITAS_USE_MPI OFF)
set(CELERITAS_USE_SWIG OFF)

FetchContent_Declare(
  adept
  GIT_REPOSITORY https://github.com/apt-sim/adept.git
  GIT_TAG        8e7b9b49f1730f8db6ff3520942fd5ae64760b21 # master, so will need to bring in code from MR 268!
)
FetchContent_Declare(
  celeritas
  GIT_REPOSITORY https://github.com/celeritas-project/celeritas.git
  GIT_TAG        c0f4129c391a6d7cfad41f3e0e3ccd0445633732 # v0.4.1
)
FetchContent_MakeAvailable(adept celeritas)

# Post fetchcontent setup
list(APPEND CMAKE_MODULE_PATH ${celeritas_SOURCE_DIR}/cmake)
find_package(Geant4 REQUIRED)

# -- Build the offloader library
include(CeleritasLibrary)
celeritas_add_library(GPUOffload
  GPUOffload/GPUOffload.cc
  GPUOffload/celeritas/CeleritasOffload.cc
  GPUOffload/celeritas/Celeritas.cc)
celeritas_target_include_directories(GPUOffload
  PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/GPUOffload>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
celeritas_target_link_libraries(GPUOffload
  Celeritas::accel
  ${Geant4_LIBRARIES}
)

# Minimal exe. Use celeritas linking to resolve RDC linking if shared libs are used
add_executable(trackingmanager-offload trackingmanager-offload.cc)
celeritas_target_link_libraries(trackingmanager-offload
  Celeritas::accel
  GPUOffload
  ${Geant4_LIBRARIES}
)
