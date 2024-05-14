#include "CeleritasOffload.hh"

#include <accel/TrackingManagerOffload.hh>

#include "Celeritas.hh"

void CeleritasOffload::Build()
{
    CelerSimpleOffload().Build(
        &CelerSetupOptions(), &CelerSharedParams(), &CelerLocalTransporter());
}

void CeleritasOffload::BuildForMaster()
{
    CelerSimpleOffload().BuildForMaster(&CelerSetupOptions(),
                                        &CelerSharedParams());
}

std::unique_ptr<G4VTrackingManager> CeleritasOffload::MakeTrackingManager()
{
    return std::make_unique<celeritas::TrackingManagerOffload>(
        &CelerSharedParams(), &CelerLocalTransporter());
}

void CeleritasOffload::BeginOfEventAction(G4Event const* event)
{
    CelerSimpleOffload().BeginOfEventAction(event);
}

void CeleritasOffload::EndOfEventAction(G4Event const* event)
{
    CelerSimpleOffload().EndOfEventAction(event);
}

void CeleritasOffload::BeginOfRunAction(G4Run const* event)
{
    CelerSimpleOffload().BeginOfRunAction(event);
}

void CeleritasOffload::EndOfRunAction(G4Run const* event)
{
    CelerSimpleOffload().EndOfRunAction(event);
}
