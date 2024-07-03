#include "NoOffload.hh"

#include <G4VTrackingManager.hh>

// Return concrete tracking manager for this offloader
std::unique_ptr<G4VTrackingManager> NoOffload::MakeTrackingManager()
{
    return std::unique_ptr<G4VTrackingManager>(nullptr);
}

//! Initialization of this class on a worker thread
void NoOffload::Build() {};

//! Initialization of this class on the master thread
void NoOffload::BuildForMaster() {}

// Perform any operations need on Run start
void NoOffload::BeginOfRunAction(G4Run const* /*run*/) {}

// Perform any operations needed on Event start
void NoOffload::BeginOfEventAction(G4Event const* /*event*/) {}

// Perform any operations needed on Event end;
void NoOffload::EndOfEventAction(G4Event const* /*event*/) {}

// Perform and operations needed on Run end;
void NoOffload::EndOfRunAction(G4Run const* /*run*/) {}
