cmake_minimum_required(VERSION 3.18...3.28)
# Must enable languages here so that calls after FetchContent works
project(AplusC LANGUAGES C CXX CUDA)

include(CTest)
include(FetchContent)

# -- Externals
# Here we now have an issue of possible duplication, so we
# 1. Setup Celeritas first,
# 2. then G4VG so it uses local Celeritas::geocel,
# 3. then AdePT so it uses local G4VG
# If/when G4VG becomes standalone, this will continue to work depending
# on how AdePT/Celeritas pick it up

# We have to set this so that packages found/required by Celeritas
# are visible when linking Celeritas to g4vg/AdePT
set(CMAKE_FIND_PACKAGE_TARGETS_GLOBAL ON)

# Set up default options here (TBD if we want users to set these)
# Or.. in CMake setting files
set(CELERITAS_USE_CUDA ON)
set(CELERITAS_USE_Geant4 ON)
set(CELERITAS_USE_VecGeom ON)
set(CELERITAS_USE_MPI OFF)
set(CELERITAS_USE_SWIG OFF)
set(CELERITAS_USE_Python OFF)
set(CELERITAS_USE_ROOT OFF)
# These are required for correct interworking...
set(CELERITAS_CORE_GEO VecGeom)
set(CELERITAS_UNITS CLHEP)
set(CELERITAS_REAL_TYPE double)

FetchContent_Declare(
  celeritas
  GIT_REPOSITORY https://github.com/celeritas-project/celeritas.git
  GIT_TAG        17d96db867322102a7bb6f18a46ecda13871cd37 # develop
)
FetchContent_MakeAvailable(celeritas)

# Force Celeritas "found" so g4vg will use Celeritas::geocel from local FetchContent
set(Celeritas_FOUND ON)
FetchContent_Declare(
  g4vg
  GIT_REPOSITORY https://github.com/celeritas-project/g4vg.git
  GIT_TAG        9fd87968a8248e0050e831fffdd89e0346b96147 # main
)
FetchContent_MakeAvailable(g4vg)

# As above, force G4VG "found" so AdePT will use G4VG from local FetchContent
# TODO: work out why this fails with "No known features for CUDA compiler" for
# the flush try-compile on first run, but succeeds on second.
set(G4VG_FOUND ON)
FetchContent_Declare(
  adept
  GIT_REPOSITORY https://github.com/drbenmorgan/adept.git
  GIT_TAG        c30469b0f074f9015d809817363c6a44c665d890 # update-g4vg
)
FetchContent_MakeAvailable(adept)

# Post fetchcontent setup
list(APPEND CMAKE_MODULE_PATH ${celeritas_SOURCE_DIR}/cmake)
find_package(Geant4 REQUIRED)

# -- Build the offloader library
include(CudaRdcUtils)
cuda_rdc_add_library(GPUOffload
  GPUOffload/GPUOffload.cc
  GPUOffload/adept/AdeptOffload.cc
  GPUOffload/celeritas/CeleritasOffload.cc
  GPUOffload/celeritas/CeleritasOffloadBuilder.cc
  GPUOffload/celeritas/Celeritas.cc
  GPUOffload/none/NoOffload.cc)
cuda_rdc_target_include_directories(GPUOffload
  PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/GPUOffload>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
cuda_rdc_target_link_libraries(GPUOffload
  Celeritas::accel
  AdePT::AdePT_G4_integration
  ${Geant4_LIBRARIES}
)

# Minimal exe. Use celeritas linking to resolve RDC linking if shared libs are used
add_executable(trackingmanager-offload trackingmanager-offload.cc)
cuda_rdc_target_link_libraries(trackingmanager-offload
  GPUOffload
  ${Geant4_LIBRARIES}
)
add_test(NAME trackingmanager-offload COMMAND trackingmanager-offload)
