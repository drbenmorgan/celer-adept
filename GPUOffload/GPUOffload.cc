#include "GPUOffload.hh"

#include <memory>
#include <G4Exception.hh>
#include <G4Threading.hh>

#include "adept/AdeptOffload.hh"
#include "celeritas/CeleritasOffloadBuilder.hh"
#include "none/NoOffload.hh"

namespace
{
// Global options for the offloader backend
GPUOffloadOptions& global_offload_options()
{
    static GPUOffloadOptions opts;
    return opts;
}

// Read only access to offloader backend options
GPUOffloadOptions const& offload_options()
{
    return global_offload_options();
}

// Build a concrete offloader from options
std::unique_ptr<GPUOffloadInterface>
build_offloader(GPUOffloadOptions const& op)
{
    // Input options must be valid
    if (!op)
    {
        G4Exception("GPUOffload.cc::build_offloader",
                    "GPUOffload0003",
                    FatalException,
                    "attempted to build GPUOffloadInterface from invalid "
                    "options");
    }

    switch (op.backend)
    {
        // Should GPUOffloadOptions need to be passed to the concrete
        // offloader or builder function, the returns below can be from calls
        // dedicated functions
        case GPUOffloadBackend::Adept:
            return std::make_unique<AdeptOffload>();
        case GPUOffloadBackend::Celeritas:
        case GPUOffloadBackend::CeleritasCPU:
            return BuildCeleritasOffload(op); 
        case GPUOffloadBackend::None:
            return std::make_unique<NoOffload>();
        default:
            // shouldn't get here...
            break;
    }
    // only logical return here, but we should never get here
    return nullptr;
}

// Return ref to thread-local offloader backend
GPUOffloadInterface& local_offloader()
{
    static G4ThreadLocal auto op = build_offloader(offload_options());
    // This assumes that op will hold be nullptr...
    return *(op.get());
}
}  // namespace

// Setup offloader backend backed. Can only be called once per process
void GPUOffloadSetup(GPUOffloadOptions const& opts)
{
    static std::mutex m;
    std::lock_guard<std::mutex> scoped_lock{m};
    GPUOffloadOptions& go = global_offload_options();

    if (go)
    {
        // Could also throw std::exception and let client handle...
        G4Exception("GPUOffloadSetup",
                    "GPUOffload0001",
                    JustWarning,
                    "GPUOffloadSetup already called");
        return;
    }

    // Input options must be valid
    if (!opts)
    {
        G4Exception("GPUOffloadSetup",
                    "GPUOffload0002",
                    FatalException,
                    "invalid GPUOffloadOptions instance passed");
    }

    go = opts;
}

GPUOffloadInterface& GPUOffload()
{
    return local_offloader();
}
