#include "GPUOffload.hh"

#include "celeritas/CeleritasOffload.hh"

#include "G4Threading.hh"
#include "adept/AdeptOffload.hh"
#include "none/NoOffload.hh"

GPUOffloadInterface& GPUOffload()
{
    // TODO: How to make the choice configurable at runtime given
    // thread locality, type, and that "GPUOffload()" is a global
    // function called in several places (so can't just use argument)
    static G4ThreadLocal CeleritasOffload co;
    // static G4ThreadLocal AdeptOffload co;
    // static G4ThreadLocal NoOffload co;
    return co;
}
