#include "GPUOffload.hh"

#include "celeritas/CeleritasOffload.hh"

#include "G4Threading.hh"

GPUOffloadInterface& GPUOffload()
{
  static G4ThreadLocal CeleritasOffload co;
  return co;
}
