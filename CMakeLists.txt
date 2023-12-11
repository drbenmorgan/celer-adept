cmake_minimum_required(VERSION 3.18...3.28)
project(AplusC)

include(FetchContent)

# Try both as FetchContent. Probably AdePT first as it's more restrictive?
# Also need G4HepEM as submodule for AdePT
# ... Add hard options here (e.g. CELERITAS_USE_Geant4)
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

