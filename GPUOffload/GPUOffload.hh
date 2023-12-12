#pragma once

#include "GPUOffloadInterface.hh"

// Offloader instance is thread local (for Celeritas only at present)
// NB: Choosing between AdePT/Celeritas might be tricky because it's required
// that the returned reference is thread-local.
GPUOffloadInterface& GPUOffload(/* argument to choose implementation ?*/);
