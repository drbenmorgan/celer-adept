#pragma once

#include "GPUOffloadInterface.hh"

class CeleritasOffload final : public GPUOffloadInterface
{
  public:
    // Return concrete tracking manager for this offloader
    std::unique_ptr<G4VTrackingManager> MakeTrackingManager() override;

    //! Initialization of this class on a worker thread
    void Build() override;

    //! Initialization of this class on the master thread
    void BuildForMaster() override;

    // Perform any operations need on Run start
    void BeginOfRunAction(G4Run const* run) override;

    // Perform any operations needed on Event start
    void BeginOfEventAction(G4Event const* event) override;

    // Perform any operations needed on Event end;
    void EndOfEventAction(G4Event const* event) override;

    // Perform and operations needed on Run end;
    void EndOfRunAction(G4Run const* run) override;
};
