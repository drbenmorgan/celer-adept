#pragma once

#include <memory>
#include <G4VTrackingManager.hh>

class G4Event;
class G4Run;

// This largely reproduces Celeritas' SimpleOffload interface, making the
// calls pure virtual so we can override these as needed for Celeritas AND
// AdePT. It adds the MakeTrackingManager to return the concrete
// G4VTrackingManager for the offload.
//
// TODO: There will certainly be setup options to pass, but this is to be
// determined what's common/divergent between AdePT/Celeritas. We try to hide
// the implementation of these AFAP.
class GPUOffloadInterface
{
  public:
    virtual ~GPUOffloadInterface() = default;

    // Return concrete tracking manager for this offloader
    virtual std::unique_ptr<G4VTrackingManager> MakeTrackingManager() = 0;

    //! Initialization of this class on a worker thread
    virtual void Build() = 0;

    //! Initialization of this class on the master thread
    virtual void BuildForMaster() = 0;

    // Perform any operations need on Run start
    virtual void BeginOfRunAction(G4Run const* run) = 0;

    // Perform any operations needed on Event start
    virtual void BeginOfEventAction(G4Event const* event) = 0;

    // Perform any operations needed on Event end;
    virtual void EndOfEventAction(G4Event const* event) = 0;

    // Perform and operations needed on Run end;
    virtual void EndOfRunAction(G4Run const* run) = 0;
};
