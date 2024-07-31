#include "AdeptCPUOffload.hh"

#include "AdeptCPUTrackingManager.hh"

//! Return concrete tracking manager for this offloader
std::unique_ptr<G4VTrackingManager> AdeptCPUOffload::MakeTrackingManager()
{
    // Modelled after AdePT's AdePTPhysics::ConstructProcess(), but adapted
    // as we don't store the transporter/config in the physics list.
    auto man = std::make_unique<AdeptCPUTrackingManager>();
    return man;
}

//! Initialization of this class on a worker thread
void AdeptCPUOffload::Build() {}

//! Initialization of this class on the master thread
void AdeptCPUOffload::BuildForMaster() {}

// Perform any operations need on Run start
void AdeptCPUOffload::BeginOfRunAction(G4Run const* /*run*/) {}

// Perform any operations needed on Event start
void AdeptCPUOffload::BeginOfEventAction(G4Event const* /*event*/) {}

// Perform any operations needed on Event end;
void AdeptCPUOffload::EndOfEventAction(G4Event const* /*event*/) {}

// Perform and operations needed on Run end;
void AdeptCPUOffload::EndOfRunAction(G4Run const* /*run*/) {}
