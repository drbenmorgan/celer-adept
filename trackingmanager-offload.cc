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
#include <G4EmParameters.hh>
#include <G4EmStandardPhysics.hh>
#include <G4Gamma.hh>
#include <G4HadronicProcessStore.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4MultiEventAction.hh>
#include <G4MultiFunctionalDetector.hh>
#include <G4MultiRunAction.hh>
#include <G4PSEnergyDeposit.hh>
#include <G4PVPlacement.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4Positron.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4Run.hh>
#include <G4RunManagerFactory.hh>
#include <G4SDManager.hh>
#include <G4StatDouble.hh>
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
class GPUOffloadRunAction final : public G4UserRunAction
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
// - Run used to process the hits/mfds
class MFDRun final : public G4Run
{
  public:
    void RecordEvent(G4Event const* event) override
    {
        auto* pHCE = event->GetHCofThisEvent();
        if (pHCE == nullptr)
            return;

        auto* pSDMan = G4SDManager::GetSDMpointer();
        auto edepID = pSDMan->GetCollectionID("scorers/edep");
        auto edepHitsMap
            = dynamic_cast<G4THitsMap<G4double>*>(pHCE->GetHC(edepID));
        if (edepHitsMap != nullptr)
        {
            for (auto& mapElement : (*edepHitsMap->GetMap()))
            {
                energy_deposition += *(mapElement.second);
            }
        }
    }

    void Merge(G4Run const* run) override
    {
        auto in = dynamic_cast<MFDRun const*>(run);
        energy_deposition += in->energy_deposition;
        G4Run::Merge(run);
    }

    G4StatDouble const& GetEnergyDeposition() const
    {
        return energy_deposition;
    }

  private:
    G4StatDouble energy_deposition;
};

// - RunAction used to process the hits/mfds
class MFDProcessingRunAction final : public G4UserRunAction
{
  public:
    void EndOfRunAction(G4Run const* run) final
    {
        if (IsMaster())
        {
            auto in = dynamic_cast<MFDRun const*>(run);
            auto average_edep = in->GetEnergyDeposition();
            G4cout << "Average energy deposition per event = "
                   << average_edep.mean() / GeV << " pm "
                   << average_edep.rms() / GeV << G4endl;
        }
    }

    G4Run* GenerateRun() { return new MFDRun; }
};

//---------------------------------------------------------------------------//
// - Combined RunAction
class RunAction final : public G4MultiRunAction
{
  public:
    RunAction()
    {
        this->emplace_back(new GPUOffloadRunAction);
        this->emplace_back(new MFDProcessingRunAction);
    }
};

//---------------------------------------------------------------------------//
// - EventAction used to notify the offloader of the Begin/End of Event
class GPUOffloadEventAction final : public G4UserEventAction
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
// - Combined EventAction
class EventAction final : public G4MultiEventAction
{
  public:
    EventAction() { this->emplace_back(new GPUOffloadEventAction); }
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
        auto mfd = new G4MultiFunctionalDetector("scorers");
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

    void GeneratePrimaries(G4Event* event) final
    {
        gun_.GeneratePrimaryVertex(event);
    }

  private:
    G4ParticleGun gun_;
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

int main()
{
    std::unique_ptr<G4RunManager> run_manager{
        G4RunManagerFactory::CreateRunManager()};
    run_manager->SetUserInitialization(new DetectorConstruction{});

    // Use FTFP_BERT, but replace EM constructor with our own that
    // overrides ConstructProcess to use GPU tracking for e-/e+/g
    // Could also be dedicated PhysList, but this is easiest demo for now.
    auto physics_list = new FTFP_BERT{/* verbosity = */ 0};
    physics_list->ReplacePhysics(new EMPhysicsConstructor);
    run_manager->SetUserInitialization(physics_list);

    // Quieten down physics so we can see what's going on
    G4EmParameters::Instance()->SetVerbose(0);
    G4HadronicProcessStore::Instance()->SetVerbose(0);

    run_manager->SetUserInitialization(new ActionInitialization());
    run_manager->Initialize();
    // Do two runs to check behaviour between
    run_manager->BeamOn(10);
    // This causes an exception in Celeritas "celeritas::activate_device may be
    // called only once per application" AdePT is fine
    //run_manager->BeamOn(16);

    return 0;
}
