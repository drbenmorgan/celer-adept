#include "AdeptOffload.hh"


// Look at AdePT's new Example22 for how to implement the functions.
// Just mocked in for now.
#include "G4VTrackingManager.hh"

// Return concrete tracking manager for this offloader
std::unique_ptr<G4VTrackingManager> AdeptOffload::MakeTrackingManager()
{
  return std::unique_ptr<G4VTrackingManager>();
}

//! Initialization of this class on a worker thread
void AdeptOffload::Build()
{};

//! Initialization of this class on the master thread
void AdeptOffload::BuildForMaster()
{}

// Perform any operations need on Run start
void AdeptOffload::BeginOfRunAction(G4Run const* /*run*/)
{}

// Perform any operations needed on Event start
void AdeptOffload::BeginOfEventAction(G4Event const* /*event*/)
{}

// Perform any operations needed on Event end;
void AdeptOffload::EndOfEventAction(G4Event const* /*event*/)
{}

// Perform and operations needed on Run end;
void AdeptOffload::EndOfRunAction(G4Run const* /*run*/)
{}
