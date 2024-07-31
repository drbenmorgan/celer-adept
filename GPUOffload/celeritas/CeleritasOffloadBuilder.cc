#include "CeleritasOffloadBuilder.hh"

#include <G4Exception.hh>
#include <corecel/Assert.hh>
#include <corecel/Macros.hh>
#include <corecel/sys/Environment.hh>

std::unique_ptr<CeleritasOffload>
BuildCeleritasOffload(GPUOffloadOptions const& opts)
{
    if (opts.backend == GPUOffloadBackend::Celeritas)
    {
        if (!CELER_USE_DEVICE)
        {
            G4Exception("CeleritasOffloadBuilder.cc::BuildCeleritasOffload",
                        "CeleritasOffloadBuilder0001",
                        FatalException,
                        "Requested Celeritas GPU offload, but Celeritas was "
                        "not compiled with GPU support");
        }
        if (!celeritas::getenv("CELER_DISABLE_DEVICE").empty())
        {
            G4Exception("CeleritasOffloadBuilder.cc::BuildCeleritasOffload",
                        "CeleritasOffloadBuilder0002",
                        FatalException,
                        "Requested Celeritas GPU offload, but "
                        "CELER_DISABLE_DEVICE environment variable set");
        }
    }
    else if (opts.backend == GPUOffloadBackend::CeleritasCPU)
    {
        // Temporary hack. Follow https://github.com/celeritas-project/celeritas/issues/1037
        // for future solution through parameters
        celeritas::environment().insert({"CELER_DISABLE_DEVICE", "1"});
        CELER_ENSURE(celeritas::getenv("CELER_DISABLE_DEVICE") == "1");
    }
    else
    {
        G4Exception("CeleritasOffloadBuilder.cc::BuildCeleritasOffload",
                    "CeleritasOffloadBuilder0003",
                    FatalException,
                    "GPUOffloadOptions.backend was not Celeritas or "
                    "CeleritasCPU");
    }
    return std::make_unique<CeleritasOffload>();
}
