#include "AdeptOffload.hh"

#include <AdePT/core/AdePTConfiguration.hh>
#include <AdePT/core/AdePTTransport.h>
#include <AdePT/integration/AdePTGeant4Integration.hh>
#include <AdePT/integration/AdePTTrackingManager.hh>

// Thread local transporter
AdePTTransport<AdePTGeant4Integration>* GetTransporter()
{
    // TODO: If created on stack/as unique_ptr (so deleted by tls at exit) get
    // segfaults..
    static G4ThreadLocal auto* transporter
        = new AdePTTransport<AdePTGeant4Integration>;
    return transporter;
}

// Global config
AdePTConfiguration& GetConfiguration()
{
    static AdePTConfiguration config;
    return config;
}

// Return concrete tracking manager for this offloader
std::unique_ptr<G4VTrackingManager> AdeptOffload::MakeTrackingManager()
{
    // Modelled after AdePT's AdePTPhysics::ConstructProcess(), but adapted
    // as we don't store the transporter/config in the physics list.
    auto man = std::make_unique<AdePTTrackingManager>();
    man->SetAdePTTransport(GetTransporter());
    man->SetAdePTConfiguration(&GetConfiguration());

    return man;
}

//! Initialization of this class on a worker thread
void AdeptOffload::Build() {}

//! Initialization of this class on the master thread
void AdeptOffload::BuildForMaster() {}

// Perform any operations need on Run start
void AdeptOffload::BeginOfRunAction(G4Run const* /*run*/) {}

// Perform any operations needed on Event start
void AdeptOffload::BeginOfEventAction(G4Event const* /*event*/) {}

// Perform any operations needed on Event end;
void AdeptOffload::EndOfEventAction(G4Event const* /*event*/) {}

// Perform and operations needed on Run end;
void AdeptOffload::EndOfRunAction(G4Run const* /*run*/) {}
