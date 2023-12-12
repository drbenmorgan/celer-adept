#pragma once

#include "GPUOffloadInterface.hh"

// Offloader instance is thread local (for Celeritas only at present)
// NB: Choosing between AdePT/Celeritas might be tricky because we
// do have a static thread_local under the hood. 
GPUOffloadInterface& GPUOffload(/* argument to choose implementation ?*/);
