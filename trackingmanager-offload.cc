//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file trackingmanager-offload.cc
//---------------------------------------------------------------------------//

#include <algorithm>
#include <iterator>
#include <type_traits>
#include <FTFP_BERT.hh>
#include <G4Box.hh>
#include <G4Electron.hh>
#include <G4EmStandardPhysics.hh>
#include <G4Gamma.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4MultiFunctionalDetector.hh>
#include <G4PSEnergyDeposit.hh>
#include <G4PVPlacement.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4Positron.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4RunManagerFactory.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Threading.hh>
#include <G4ThreeVector.hh>
#include <G4Track.hh>
#include <G4TrackStatus.hh>
#include <G4Types.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserActionInitialization.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>


//---------------------------------------------------------------------------//
// - Parts coupled to the offload

// Frontend "factory" function
#include "GPUOffload/GPUOffload.hh"

//---------------------------------------------------------------------------//
// - RunAction used to notify the offloader of the Begin/End of Run
class RunAction final : public G4UserRunAction
{
  public:
    void BeginOfRunAction(G4Run const* run) final
    {
        GPUOffload().BeginOfRunAction(run);
    }
    void EndOfRunAction(G4Run const* run) final
    {
        GPUOffload().EndOfRunAction(run);
    }
};

//---------------------------------------------------------------------------//
// - EventAction used to notify the offloader of the Begin/End of Event
class EventAction final : public G4UserEventAction
{
  public:
    void BeginOfEventAction(G4Event const* event) final
    {
        GPUOffload().BeginOfEventAction(event);
    }

    void EndOfEventAction(G4Event const* event) final
    {
        GPUOffload().EndOfEventAction(event);
    }
};

//---------------------------------------------------------------------------//
// - Set up the needed Run/Event Actions
// - Notify the offloader of the Build/BuildForMaster states
class ActionInitialization final : public G4VUserActionInitialization
{
  public:
    void BuildForMaster() const final
    {
        GPUOffload().BuildForMaster();
        this->SetUserAction(new RunAction{});
    }
    void Build() const final
    {
        GPUOffload().Build();
        this->SetUserAction(new PrimaryGeneratorAction{});
        this->SetUserAction(new RunAction{});
        this->SetUserAction(new EventAction{});
    }
};


//---------------------------------------------------------------------------//
// - Custom PhysicsConstructor to add the Offloader's G4VTrackingManager to
//   e-/e+/gamma particles.
// - Inherited from G4EmStandardPhysics because we want to ensure that the
//   CPU processes/models are constructed properly (To be reviewed, and see
//   also G4 RE07 example).
class EMPhysicsConstructor final : public G4EmStandardPhysics
{
  public:
    using G4EmStandardPhysics::G4EmStandardPhysics;

    void ConstructProcess() override
    {
        G4EmStandardPhysics::ConstructProcess();

        // Create and add the chosen GPU tracking manager to EM particles
        // NB: this ultimately should be thread-local, but G4 -should- ensure
        // that for us
        auto gpu_tracking = GPUOffload().MakeTrackingManager();
        G4Electron::Definition()->SetTrackingManager(gpu_tracking.get());
        G4Positron::Definition()->SetTrackingManager(gpu_tracking.get());
        G4Gamma::Definition()->SetTrackingManager(gpu_tracking.get());
        // Geant4 will take ownership of the pointer, so we can release it here
        gpu_tracking.release();
    }
};

//---------------------------------------------------------------------------//
// Non-coupled to offloader as yet
//---------------------------------------------------------------------------//
class DetectorConstruction final : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction()
        : aluminum_{new G4Material{
            "Aluminium", 13., 26.98 * g / mole, 2.700 * g / cm3}}
    {
    }

    G4VPhysicalVolume* Construct() final
    {
        auto* box = new G4Box("world", 1000 * cm, 1000 * cm, 1000 * cm);
        auto* lv = new G4LogicalVolume(box, aluminum_, "world");
        auto* pv = new G4PVPlacement(
            0, G4ThreeVector{}, lv, "world", nullptr, false, 0);
        return pv;
    }

    // Celeritas will fail if sds are active but none are found, so we just
    // add a simple PS here (we don't use the results yet)
    void ConstructSDandField()
    {
        auto mfd = new G4MultiFunctionalDetector("world");
        G4SDManager::GetSDMpointer()->AddNewDetector(mfd);
        auto scorer = new G4PSEnergyDeposit("edep");
        mfd->RegisterPrimitive(scorer);
        SetSensitiveDetector("world", mfd);
    }

  private:
    G4Material* aluminum_;
};

//---------------------------------------------------------------------------//
class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction()
    {
        auto g4particle_def
            = G4ParticleTable::GetParticleTable()->FindParticle(2112);
        gun_.SetParticleDefinition(g4particle_def);
        gun_.SetParticleEnergy(100 * GeV);
        gun_.SetParticlePosition(G4ThreeVector{0, 0, 0});  // origin
        gun_.SetParticleMomentumDirection(G4ThreeVector{1, 0, 0});  // +x
    }

    // Generate 100 GeV neutrons
    void GeneratePrimaries(G4Event* event) final
    {
        gun_.GeneratePrimaryVertex(event);
    }

  private:
    G4ParticleGun gun_;
};


//---------------------------------------------------------------------------//

int main()
{
    std::unique_ptr<G4RunManager> run_manager{
        G4RunManagerFactory::CreateRunManager()};  // G4RunManagerType::SerialOnly)};

    run_manager->SetUserInitialization(new DetectorConstruction{});

    // Use FTFP_BERT, but replace EM constructor with our own that
    // overrides ConstructProcess to use GPU tracking for e-/e+/g
    // Could also be dedicated PhysList, but this is easiest demo for now.
    auto physics_list = new FTFP_BERT{/* verbosity = */ 0};
    physics_list->ReplacePhysics(new EMPhysicsConstructor);
    run_manager->SetUserInitialization(physics_list);

    run_manager->SetUserInitialization(new ActionInitialization());
    run_manager->Initialize();
    // Do two runs to check behaviour between
    run_manager->BeamOn(8);
    // This causes an exception in Celeritas "celeritas::activate_device may be called only once per application"
    //run_manager->BeamOn(16);

    return 0;
}
