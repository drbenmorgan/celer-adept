#include "AdeptCPUTrackingManager.hh"

#include <G4Electron.hh>
#include <G4EmParameters.hh>
#include <G4EmProcessSubType.hh>
#include <G4Gamma.hh>
#include <G4HepEmData.hh>
#include <G4HepEmElectronManager.hh>
#include <G4HepEmElectronTrack.hh>
#include <G4HepEmGammaManager.hh>
#include <G4HepEmGammaTrack.hh>
#include <G4HepEmMatCutData.hh>
#include <G4HepEmNoProcess.hh>
#include <G4HepEmPositronInteractionAnnihilation.hh>
#include <G4HepEmRandomEngine.hh>
#include <G4HepEmRunManager.hh>
#include <G4HepEmTLData.hh>
#include <G4MaterialCutsCouple.hh>
#include <G4Positron.hh>
#include <G4ProcessType.hh>
#include <G4ProductionCutsTable.hh>
#include <G4SafetyHelper.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4StepStatus.hh>
#include <G4Threading.hh>
#include <G4Track.hh>
#include <G4TransportationManager.hh>
#include <G4TransportationProcessType.hh>

#include "TrackingManagerHelper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AdeptCPUTrackingManager::AdeptCPUTrackingManager()
{
    fRunManager = new G4HepEmRunManager(G4Threading::IsMasterThread());
    fRandomEngine = new G4HepEmRandomEngine(G4Random::getTheEngine());
    fSafetyHelper
        = G4TransportationManager::GetTransportationManager()->GetSafetyHelper();
    fSafetyHelper->InitialiseHelper();
    fStep = new G4Step;
    fStep->NewSecondaryVector();

    // Construct fake G4VProcess-es with the proper name and indices matching
    // the hepEm process indices
    fElectronNoProcessVector.push_back(
        new G4HepEmNoProcess("eIoni",
                             G4ProcessType::fElectromagnetic,
                             G4EmProcessSubType::fIonisation));
    fElectronNoProcessVector.push_back(
        new G4HepEmNoProcess("eBrem",
                             G4ProcessType::fElectromagnetic,
                             G4EmProcessSubType::fBremsstrahlung));
    fElectronNoProcessVector.push_back(
        new G4HepEmNoProcess("annihl",
                             G4ProcessType::fElectromagnetic,
                             G4EmProcessSubType::fAnnihilation));
    fElectronNoProcessVector.push_back(
        new G4HepEmNoProcess("msc",
                             G4ProcessType::fElectromagnetic,
                             G4EmProcessSubType::fMultipleScattering));
    fGammaNoProcessVector.push_back(
        new G4HepEmNoProcess("conv",
                             G4ProcessType::fElectromagnetic,
                             G4EmProcessSubType::fGammaConversion));
    fGammaNoProcessVector.push_back(
        new G4HepEmNoProcess("compt",
                             G4ProcessType::fElectromagnetic,
                             G4EmProcessSubType::fComptonScattering));
    fGammaNoProcessVector.push_back(
        new G4HepEmNoProcess("phot",
                             G4ProcessType::fElectromagnetic,
                             G4EmProcessSubType::fPhotoElectricEffect));
    fTransportNoProcess = new G4HepEmNoProcess(
        "Transportation", G4ProcessType::fTransportation, TRANSPORTATION);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AdeptCPUTrackingManager::~AdeptCPUTrackingManager()
{
    // Per behaviour in Physics Constructors, we do not delete the
    // G4HepEmNoProcess instances as these are owned by G4ProcessTable.
    delete fRunManager;
    delete fRandomEngine;
    delete fStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AdeptCPUTrackingManager::BuildPhysicsTable(G4ParticleDefinition const& part)
{
    if (&part == G4Electron::Definition())
    {
        fRunManager->Initialize(fRandomEngine, 0);
    }
    else if (&part == G4Positron::Definition())
    {
        fRunManager->Initialize(fRandomEngine, 1);
    }
    else if (&part == G4Gamma::Definition())
    {
        fRunManager->Initialize(fRandomEngine, 2);
    }
    else
    {
        std::cerr << " **** ERROR in G4HepEmProcess::BuildPhysicsTable: "
                     "unknown particle "
                  << std::endl;
        exit(-1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AdeptCPUTrackingManager::PreparePhysicsTable(
    G4ParticleDefinition const& part)
{
    applyCuts = G4EmParameters::Instance()->ApplyCuts();

    if (applyCuts)
    {
        auto* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
        theCutsGamma = theCoupleTable->GetEnergyCutsVector(idxG4GammaCut);
        theCutsElectron = theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut);
        theCutsPositron = theCoupleTable->GetEnergyCutsVector(idxG4PositronCut);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AdeptCPUTrackingManager::TrackElectron(G4Track* aTrack)
{
    TrackingManagerHelper::ChargedNavigation navigation;

    // Prepare for calling the user action.
    auto* evtMgr = G4EventManager::GetEventManager();
    auto* userTrackingAction = evtMgr->GetUserTrackingAction();
    auto* userSteppingAction = evtMgr->GetUserSteppingAction();

    // Locate the track in geometry.
    {
        auto* transMgr = G4TransportationManager::GetTransportationManager();
        auto* linearNavigator = transMgr->GetNavigatorForTracking();

        G4ThreeVector const& pos = aTrack->GetPosition();
        G4ThreeVector const& dir = aTrack->GetMomentumDirection();

        // Do not assign directly, doesn't work if the handle is empty.
        G4TouchableHandle touchableHandle;
        if (aTrack->GetTouchableHandle())
        {
            touchableHandle = aTrack->GetTouchableHandle();
            // FIXME: This assumes we only ever have G4TouchableHistorys!
            auto* touchableHistory = (G4TouchableHistory*)touchableHandle();
            G4VPhysicalVolume* oldTopVolume = touchableHandle->GetVolume();
            G4VPhysicalVolume* newTopVolume
                = linearNavigator->ResetHierarchyAndLocate(
                    pos, dir, *touchableHistory);
            // TODO: WHY?!
            if (newTopVolume != oldTopVolume
                || oldTopVolume->GetRegularStructureId() == 1)
            {
                touchableHandle = linearNavigator->CreateTouchableHistory();
                aTrack->SetTouchableHandle(touchableHandle);
            }
        }
        else
        {
            linearNavigator->LocateGlobalPointAndSetup(pos, &dir, false, false);
            touchableHandle = linearNavigator->CreateTouchableHistory();
            aTrack->SetTouchableHandle(touchableHandle);
        }
        aTrack->SetNextTouchableHandle(touchableHandle);
    }

    // Prepare data structures used while tracking.
    G4Step& step = *fStep;
    G4TrackVector& secondaries = *step.GetfSecondary();
    G4StepPoint& preStepPoint = *step.GetPreStepPoint();
    G4StepPoint& postStepPoint = *step.GetPostStepPoint();
    step.InitializeStep(aTrack);
    aTrack->SetStep(&step);

    // Start of tracking: Inform user and processes.
    if (userTrackingAction)
    {
        userTrackingAction->PreUserTrackingAction(aTrack);
    }

    // === StartTracking ===
    G4HepEmTLData* theTLData = fRunManager->GetTheTLData();
    G4HepEmElectronTrack* theElTrack = theTLData->GetPrimaryElectronTrack();
    G4HepEmTrack* thePrimaryTrack = theElTrack->GetTrack();
    theElTrack->ReSet();
    // In principle, we could continue to use the other generated Gaussian
    // number as long as we are in the same event, but play it safe.
    G4HepEmRandomEngine* rnge = theTLData->GetRNGEngine();
    rnge->DiscardGauss();

    // Pull data structures into local variables.
    G4HepEmData* theHepEmData = fRunManager->GetHepEmData();
    G4HepEmParameters* theHepEmPars = fRunManager->GetHepEmParameters();

    G4DynamicParticle const* theG4DPart = aTrack->GetDynamicParticle();
    G4int const trackID = aTrack->GetTrackID();

    // Init state that never changes for a track.
    double const charge = aTrack->GetParticleDefinition()->GetPDGCharge();
    bool const isElectron = (charge < 0.0);
    thePrimaryTrack->SetCharge(charge);
    // === StartTracking ===

    while (aTrack->GetTrackStatus() == fAlive)
    {
        // Beginning of this step: Prepare data structures.
        aTrack->IncrementCurrentStepNumber();

        step.CopyPostToPreStepPoint();
        step.ResetTotalEnergyDeposit();
        G4TouchableHandle const& touchableHandle
            = aTrack->GetNextTouchableHandle();
        aTrack->SetTouchableHandle(touchableHandle);

        auto* lvol = aTrack->GetTouchable()->GetVolume()->GetLogicalVolume();
        preStepPoint.SetMaterial(lvol->GetMaterial());
        auto* MCC = lvol->GetMaterialCutsCouple();
        preStepPoint.SetMaterialCutsCouple(lvol->GetMaterialCutsCouple());

        // Query step lengths from pyhsics and geometry, decide on limit.
        G4double const preStepEkin = theG4DPart->GetKineticEnergy();
        G4double const preStepLogEkin = theG4DPart->GetLogKineticEnergy();
        thePrimaryTrack->SetEKin(preStepEkin, preStepLogEkin);

        int const g4IMC = MCC->GetIndex();
        int const hepEmIMC
            = theHepEmData->fTheMatCutData->fG4MCIndexToHepEmMCIndex[g4IMC];
        thePrimaryTrack->SetMCIndex(hepEmIMC);
        bool preStepOnBoundary = preStepPoint.GetStepStatus()
                                 == G4StepStatus::fGeomBoundary;
        thePrimaryTrack->SetOnBoundary(preStepOnBoundary);
        double const preSafety
            = preStepOnBoundary
                  ? 0.
                  : fSafetyHelper->ComputeSafety(aTrack->GetPosition());
        thePrimaryTrack->SetSafety(preSafety);

        // Sample the `number-of-interaction-left`
        for (int ip = 0; ip < 3; ++ip)
        {
            if (thePrimaryTrack->GetNumIALeft(ip) <= 0.)
            {
                thePrimaryTrack->SetNumIALeft(-G4HepEmLog(rnge->flat()), ip);
            }
        }
        // True distance to discrete interaction.
        G4HepEmElectronManager::HowFarToDiscreteInteraction(
            theHepEmData, theHepEmPars, theElTrack);
        // Remember which process was selected - MSC might limit the sub-steps.
        int const iDProc = thePrimaryTrack->GetWinnerProcessIndex();

        double stepLimitLeft = theElTrack->GetPStepLength();
        double totalTruePathLength = 0, totalEloss = 0;
        bool continueStepping = fMultipleSteps, stopped = false;

        theElTrack->SavePreStepEKin();

        do
        {
            // Possibly true step limit of MSC, and conversion to geometrical
            // step length.
            G4HepEmElectronManager::HowFarToMSC(
                theHepEmData, theHepEmPars, theElTrack, rnge);
            if (thePrimaryTrack->GetWinnerProcessIndex() != -2)
            {
                // If MSC did not limit the step, exit the loop after this
                // iteration.
                continueStepping = false;
            }

            // Get the geometrcal step length: straight line distance to make
            // along the original direction.
            G4double physicalStep = thePrimaryTrack->GetGStepLength();
            G4double geometryStep
                = navigation.MakeStep(*aTrack, step, physicalStep);

            bool geometryLimitedStep = geometryStep < physicalStep;
            G4double finalStep = geometryLimitedStep ? geometryStep
                                                     : physicalStep;

            step.UpdateTrack();

            navigation.FinishStep(*aTrack, step);

            if (geometryLimitedStep)
            {
                continueStepping = false;
                // Check if the track left the world.
                if (aTrack->GetNextVolume() == nullptr)
                {
                    aTrack->SetTrackStatus(fStopAndKill);
                    break;
                }
            }

            bool const postStepOnBoundary = postStepPoint.GetStepStatus()
                                            == G4StepStatus::fGeomBoundary;

            // NOTE: this primary track is the same as in the last call in the
            // HowFar()
            //       But transportation might changed its direction,
            //       geomertical step length, or status ( on boundary or not).
            G4ThreeVector const& primDir = theG4DPart->GetMomentumDirection();
            thePrimaryTrack->SetDirection(primDir[0], primDir[1], primDir[2]);
            thePrimaryTrack->SetGStepLength(finalStep);
            thePrimaryTrack->SetOnBoundary(postStepOnBoundary);
            // invoke the physics interactions (all i.e. all along- and
            // post-step as well as possible at rest)

            if (finalStep > 0)
            {
                do
                {
                    //
                    // === 1. MSC should be invoked to obtain the physics step
                    // Length
                    G4HepEmElectronManager::UpdatePStepLength(theElTrack);
                    double const pStepLength = theElTrack->GetPStepLength();
                    totalTruePathLength += pStepLength;

                    if (pStepLength <= 0.0)
                    {
                        break;
                    }
                    // compute the energy loss first based on the new step
                    // length: it will be needed in the MSC scatteirng and
                    // displacement computation here as well (that is done only
                    // if not the last step with the particle). But update the
                    // number of interaction length left before.
                    //
                    // === 2. The `number-of-interaction-left` needs to be
                    // updated based on the actual
                    //        physical step Length
                    G4HepEmElectronManager::UpdateNumIALeft(theElTrack);
                    //
                    // === 3. Continuous energy loss needs to be computed
                    stopped = G4HepEmElectronManager::ApplyMeanEnergyLoss(
                        theHepEmData, theHepEmPars, theElTrack);
                    totalEloss += thePrimaryTrack->GetEnergyDeposit();
                    if (stopped)
                    {
                        continueStepping = false;
                        break;
                    }

                    // === 4. Sample MSC direction change and displacement.
                    G4HepEmElectronManager::SampleMSC(
                        theHepEmData, theHepEmPars, theElTrack, rnge);

                    double const* pdir = thePrimaryTrack->GetDirection();
                    postStepPoint.SetMomentumDirection(
                        G4ThreeVector(pdir[0], pdir[1], pdir[2]));

                    // apply MSC displacement if its length is longer than a
                    // minimum and we are not on boundary
                    G4ThreeVector position = postStepPoint.GetPosition();
                    if (!postStepOnBoundary)
                    {
                        double const* displacement
                            = theElTrack->GetMSCTrackData()->GetDisplacement();
                        double const dLength2
                            = displacement[0] * displacement[0]
                              + displacement[1] * displacement[1]
                              + displacement[2] * displacement[2];
                        double const kGeomMinLength = 5.0e-8;  // 0.05 [nm]
                        double const kGeomMinLength2
                            = kGeomMinLength * kGeomMinLength;  // (0.05
                                                                // [nm])^2
                        if (dLength2 > kGeomMinLength2)
                        {
                            // apply displacement
                            bool isPositionChanged = true;
                            double const dispR = std::sqrt(dLength2);
                            double const postSafety
                                = 0.99
                                  * fSafetyHelper->ComputeSafety(position,
                                                                 dispR);
                            G4ThreeVector const theDisplacement(
                                displacement[0],
                                displacement[1],
                                displacement[2]);
                            // far away from geometry boundary
                            if (postSafety > 0.0 && dispR <= postSafety)
                            {
                                position += theDisplacement;
                                // near the boundary
                            }
                            else
                            {
                                // displaced point is definitely within the
                                // volume
                                if (dispR < postSafety)
                                {
                                    position += theDisplacement;
                                    // reduced displacement
                                }
                                else if (postSafety > kGeomMinLength)
                                {
                                    position += theDisplacement
                                                * (postSafety / dispR);
                                    // very small postSafety
                                }
                                else
                                {
                                    isPositionChanged = false;
                                }
                            }
                            if (isPositionChanged)
                            {
                                fSafetyHelper->ReLocateWithinVolume(position);
                                postStepPoint.SetPosition(position);
                            }
                        }
                    }

                } while (0);

                if (continueStepping)
                {
                    // Reset the energy deposit, we're accumulating it in
                    // totalEloss.
                    thePrimaryTrack->SetEnergyDeposit(0);
                    // Also reset the selected winner process that MSC
                    // replaced.
                    thePrimaryTrack->SetWinnerProcessIndex(iDProc);

                    // Save the current energy for the next MSC invocation.
                    postStepPoint.SetKineticEnergy(thePrimaryTrack->GetEKin());
                    theElTrack->SavePreStepEKin();
                    // Apply everything to the track, so that navigation sees
                    // it.
                    step.UpdateTrack();

                    // Subtract the path length we just traveled, and set this
                    // as the next attempted sub-step.
                    double const pStepLength = theElTrack->GetPStepLength();
                    stepLimitLeft -= pStepLength;
                    theElTrack->SetPStepLength(stepLimitLeft);
                    thePrimaryTrack->SetGStepLength(stepLimitLeft);

                    // Also reduce the range accordingly.
                    double const range = theElTrack->GetRange() - pStepLength;
                    theElTrack->SetRange(range);
                }
            }
        } while (continueStepping);

        // Restore the total (mean) energy loss accumulated along the
        // sub-steps.
        thePrimaryTrack->SetEnergyDeposit(totalEloss);

        if (!stopped)
        {
            // If not already stopped, restore the pre-step energy and sample
            // loss fluctuations.
            theElTrack->SetPreStepEKin(preStepEkin, preStepLogEkin);
            stopped = G4HepEmElectronManager::SampleLossFluctuations(
                theHepEmData, theHepEmPars, theElTrack, rnge);
        }

        G4VProcess const* proc = nullptr;
        if (stopped)
        {
            // call annihilation for e+ !!!
            if (!isElectron)
            {
                G4HepEmPositronInteractionAnnihilation::Perform(theTLData,
                                                                true);
                proc = fElectronNoProcessVector[2];
            }
            else
            {
                // otherwise ionization limited the step
                proc = fElectronNoProcessVector[0];
            }
        }
        else if (aTrack->GetTrackStatus() != fStopAndKill)
        {
            // === 4. Discrete part of the interaction (if any)
            G4HepEmElectronManager::PerformDiscrete(
                theHepEmData, theHepEmPars, theTLData);
            double const* pdir = thePrimaryTrack->GetDirection();
            postStepPoint.SetMomentumDirection(
                G4ThreeVector(pdir[0], pdir[1], pdir[2]));

            // Get the final process defining the step - might still be MSC!
            int const iDProc = thePrimaryTrack->GetWinnerProcessIndex();
            if (thePrimaryTrack->GetOnBoundary())
            {
                proc = fTransportNoProcess;
            }
            else if (iDProc == -1)
            {
                // ionization
                proc = fElectronNoProcessVector[0];
            }
            else if (iDProc == -2)
            {
                proc = fElectronNoProcessVector[3];
            }
            else
            {
                proc = fElectronNoProcessVector[iDProc];
            }
        }
        else
        {
            // Else the particle left the world.
            proc = fTransportNoProcess;
        }

        postStepPoint.SetProcessDefinedStep(proc);
        step.SetStepLength(totalTruePathLength);

        // energy, e-depo and status
        double const ekin = thePrimaryTrack->GetEKin();
        double edep = thePrimaryTrack->GetEnergyDeposit();
        postStepPoint.SetKineticEnergy(ekin);
        if (ekin <= 0.0)
        {
            aTrack->SetTrackStatus(fStopAndKill);
        }
        step.UpdateTrack();

        int const numSecElectron = theTLData->GetNumSecondaryElectronTrack();
        int const numSecGamma = theTLData->GetNumSecondaryGammaTrack();
        int const numSecondaries = numSecElectron + numSecGamma;
        if (numSecondaries > 0)
        {
            G4ThreeVector const& theG4PostStepPointPosition
                = postStepPoint.GetPosition();
            G4double const theG4PostStepGlobalTime
                = postStepPoint.GetGlobalTime();
            for (int is = 0; is < numSecElectron; ++is)
            {
                G4HepEmTrack* secTrack
                    = theTLData->GetSecondaryElectronTrack(is)->GetTrack();
                double const secEKin = secTrack->GetEKin();
                bool const isElectron = secTrack->GetCharge() < 0.0;
                if (applyCuts)
                {
                    if (isElectron && secEKin < (*theCutsElectron)[g4IMC])
                    {
                        edep += secEKin;
                        continue;
                    }
                    else if (!isElectron
                             && CLHEP::electron_mass_c2 < (*theCutsGamma)[g4IMC]
                             && secEKin < (*theCutsPositron)[g4IMC])
                    {
                        edep += secEKin + 2 * CLHEP::electron_mass_c2;
                        continue;
                    }
                }

                double const* dir = secTrack->GetDirection();
                G4ParticleDefinition const* partDef = G4Electron::Definition();
                if (!isElectron)
                {
                    partDef = G4Positron::Definition();
                }
                G4DynamicParticle* dp = new G4DynamicParticle(
                    partDef, G4ThreeVector(dir[0], dir[1], dir[2]), secEKin);
                G4Track* aG4Track = new G4Track(
                    dp, theG4PostStepGlobalTime, theG4PostStepPointPosition);
                aG4Track->SetParentID(trackID);
                aG4Track->SetCreatorProcess(proc);
                aG4Track->SetTouchableHandle(touchableHandle);
                secondaries.push_back(aG4Track);
            }
            theTLData->ResetNumSecondaryElectronTrack();

            for (int is = 0; is < numSecGamma; ++is)
            {
                G4HepEmTrack* secTrack
                    = theTLData->GetSecondaryGammaTrack(is)->GetTrack();
                double const secEKin = secTrack->GetEKin();
                if (applyCuts && secEKin < (*theCutsGamma)[g4IMC])
                {
                    edep += secEKin;
                    continue;
                }

                double const* dir = secTrack->GetDirection();
                G4DynamicParticle* dp = new G4DynamicParticle(
                    G4Gamma::Definition(),
                    G4ThreeVector(dir[0], dir[1], dir[2]),
                    secEKin);
                G4Track* aG4Track = new G4Track(
                    dp, theG4PostStepGlobalTime, theG4PostStepPointPosition);
                aG4Track->SetParentID(trackID);
                aG4Track->SetCreatorProcess(proc);
                aG4Track->SetTouchableHandle(touchableHandle);
                secondaries.push_back(aG4Track);
            }
            theTLData->ResetNumSecondaryGammaTrack();
        }

        step.AddTotalEnergyDeposit(edep);

        // Need to get the true step length, not the geometry step length!
        aTrack->AddTrackLength(step.GetStepLength());

        // End of this step: Call sensitive detector and stepping actions.
        if (step.GetControlFlag() != AvoidHitInvocation)
        {
            auto* sensitive = lvol->GetSensitiveDetector();
            if (sensitive)
            {
                sensitive->Hit(&step);
            }
        }

        if (userSteppingAction)
        {
            userSteppingAction->UserSteppingAction(&step);
        }

        auto* regionalAction = lvol->GetRegion()->GetRegionalSteppingAction();
        if (regionalAction)
        {
            regionalAction->UserSteppingAction(&step);
        }
    }

    // End of tracking: Inform processes and user.
    // === EndTracking ===

    if (userTrackingAction)
    {
        userTrackingAction->PostUserTrackingAction(aTrack);
    }

    evtMgr->StackTracks(&secondaries);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AdeptCPUTrackingManager::TrackGamma(G4Track* aTrack)
{
    class GammaPhysics final : public TrackingManagerHelper::Physics
    {
      public:
        GammaPhysics(AdeptCPUTrackingManager& mgr) : fMgr(mgr) {}

        void StartTracking(G4Track* aTrack) override
        {
            fMgr.fRunManager->GetTheTLData()->GetPrimaryGammaTrack()->ReSet();
            // In principle, we could continue to use the other generated
            // Gaussian number as long as we are in the same event, but play it
            // safe.
            fMgr.fRunManager->GetTheTLData()->GetRNGEngine()->DiscardGauss();
        }

        G4double GetPhysicalInteractionLength(G4Track const& track) override
        {
            G4HepEmTLData* theTLData = fMgr.fRunManager->GetTheTLData();
            G4HepEmTrack* thePrimaryTrack
                = theTLData->GetPrimaryGammaTrack()->GetTrack();
            G4HepEmData* theHepEmData = fMgr.fRunManager->GetHepEmData();
            thePrimaryTrack->SetCharge(0);
            G4DynamicParticle const* theG4DPart = track.GetDynamicParticle();
            thePrimaryTrack->SetEKin(theG4DPart->GetKineticEnergy(),
                                     theG4DPart->GetLogKineticEnergy());

            int const g4IMC = track.GetMaterialCutsCouple()->GetIndex();
            int const hepEmIMC
                = theHepEmData->fTheMatCutData->fG4MCIndexToHepEmMCIndex[g4IMC];
            thePrimaryTrack->SetMCIndex(hepEmIMC);
            G4StepPoint const* theG4PreStepPoint
                = track.GetStep()->GetPreStepPoint();
            thePrimaryTrack->SetOnBoundary(theG4PreStepPoint->GetStepStatus()
                                           == G4StepStatus::fGeomBoundary);
            G4HepEmGammaManager::HowFar(theHepEmData,
                                        fMgr.fRunManager->GetHepEmParameters(),
                                        theTLData);
            // returns with the geometrcal step length: straight line distance
            // to make along the org direction
            return thePrimaryTrack->GetGStepLength();
        }

        void
        AlongStepDoIt(G4Track& track, G4Step& step, G4TrackVector&) override
        {
            // Nothing to do here!
        }

        void PostStepDoIt(G4Track& track,
                          G4Step& step,
                          G4TrackVector& secondaries) override
        {
            G4HepEmTLData* theTLData = fMgr.fRunManager->GetTheTLData();
            G4StepPoint* theG4PostStepPoint = step.GetPostStepPoint();
            bool const onBoundary = theG4PostStepPoint->GetStepStatus()
                                    == G4StepStatus::fGeomBoundary;
            G4HepEmTrack* thePrimaryTrack
                = theTLData->GetPrimaryGammaTrack()->GetTrack();

            if (onBoundary)
            {
                thePrimaryTrack->SetGStepLength(track.GetStepLength());
                G4HepEmGammaManager::UpdateNumIALeft(thePrimaryTrack);
                theG4PostStepPoint->SetProcessDefinedStep(
                    fMgr.fTransportNoProcess);
                return;
            }
            // NOTE: this primary track is the same as in the last call in the
            // HowFar()
            //       But transportation might changed its direction,
            //       geomertical step length, or status ( on boundary or not).
            G4ThreeVector const& primDir
                = track.GetDynamicParticle()->GetMomentumDirection();
            thePrimaryTrack->SetDirection(primDir[0], primDir[1], primDir[2]);
            thePrimaryTrack->SetGStepLength(track.GetStepLength());
            thePrimaryTrack->SetOnBoundary(onBoundary);
            // invoke the physics interactions (all i.e. all along- and
            // post-step as well as possible at rest)
            G4HepEmGammaManager::Perform(fMgr.fRunManager->GetHepEmData(),
                                         fMgr.fRunManager->GetHepEmParameters(),
                                         theTLData);

            int const iDProc = thePrimaryTrack->GetWinnerProcessIndex();
            G4VProcess const* proc = fMgr.fGammaNoProcessVector[iDProc];
            theG4PostStepPoint->SetProcessDefinedStep(proc);

            // energy, e-depo, momentum direction and status
            double const ekin = thePrimaryTrack->GetEKin();
            double edep = thePrimaryTrack->GetEnergyDeposit();
            theG4PostStepPoint->SetKineticEnergy(ekin);
            if (ekin <= 0.0)
            {
                track.SetTrackStatus(fStopAndKill);
            }
            double const* pdir = thePrimaryTrack->GetDirection();
            theG4PostStepPoint->SetMomentumDirection(
                G4ThreeVector(pdir[0], pdir[1], pdir[2]));

            step.UpdateTrack();

            int const g4IMC
                = step.GetPreStepPoint()->GetMaterialCutsCouple()->GetIndex();
            // secondary: only possible is e- or gamma at the moemnt
            int const numSecElectron
                = theTLData->GetNumSecondaryElectronTrack();
            int const numSecGamma = theTLData->GetNumSecondaryGammaTrack();
            int const numSecondaries = numSecElectron + numSecGamma;
            if (numSecondaries > 0)
            {
                G4ThreeVector const& theG4PostStepPointPosition
                    = theG4PostStepPoint->GetPosition();
                G4double const theG4PostStepGlobalTime
                    = theG4PostStepPoint->GetGlobalTime();
                G4TouchableHandle const& theG4TouchableHandle
                    = track.GetTouchableHandle();
                for (int is = 0; is < numSecElectron; ++is)
                {
                    G4HepEmTrack* secTrack
                        = theTLData->GetSecondaryElectronTrack(is)->GetTrack();
                    double const secEKin = secTrack->GetEKin();
                    bool const isElectron = secTrack->GetCharge() < 0.0;
                    if (fMgr.applyCuts)
                    {
                        if (isElectron
                            && secEKin < (*fMgr.theCutsElectron)[g4IMC])
                        {
                            edep += secEKin;
                            continue;
                        }
                        else if (!isElectron
                                 && CLHEP::electron_mass_c2
                                        < (*fMgr.theCutsGamma)[g4IMC]
                                 && secEKin < (*fMgr.theCutsPositron)[g4IMC])
                        {
                            edep += secEKin + 2 * CLHEP::electron_mass_c2;
                            continue;
                        }
                    }

                    double const* dir = secTrack->GetDirection();
                    G4ParticleDefinition const* partDef
                        = G4Electron::Definition();
                    if (!isElectron)
                    {
                        partDef = G4Positron::Definition();
                    }
                    G4DynamicParticle* dp = new G4DynamicParticle(
                        partDef, G4ThreeVector(dir[0], dir[1], dir[2]), secEKin);
                    G4Track* aG4Track = new G4Track(dp,
                                                    theG4PostStepGlobalTime,
                                                    theG4PostStepPointPosition);
                    aG4Track->SetParentID(track.GetTrackID());
                    aG4Track->SetCreatorProcess(proc);
                    aG4Track->SetTouchableHandle(theG4TouchableHandle);
                    secondaries.push_back(aG4Track);
                }
                theTLData->ResetNumSecondaryElectronTrack();

                for (int is = 0; is < numSecGamma; ++is)
                {
                    G4HepEmTrack* secTrack
                        = theTLData->GetSecondaryGammaTrack(is)->GetTrack();
                    double const secEKin = secTrack->GetEKin();
                    if (fMgr.applyCuts && secEKin < (*fMgr.theCutsGamma)[g4IMC])
                    {
                        edep += secEKin;
                        continue;
                    }

                    double const* dir = secTrack->GetDirection();
                    G4DynamicParticle* dp = new G4DynamicParticle(
                        G4Gamma::Definition(),
                        G4ThreeVector(dir[0], dir[1], dir[2]),
                        secEKin);
                    G4Track* aG4Track = new G4Track(dp,
                                                    theG4PostStepGlobalTime,
                                                    theG4PostStepPointPosition);
                    aG4Track->SetParentID(track.GetTrackID());
                    aG4Track->SetCreatorProcess(proc);
                    aG4Track->SetTouchableHandle(theG4TouchableHandle);
                    secondaries.push_back(aG4Track);
                }
                theTLData->ResetNumSecondaryGammaTrack();
            }

            step.AddTotalEnergyDeposit(edep);
        }

      private:
        AdeptCPUTrackingManager& fMgr;
    };

    GammaPhysics physics(*this);
    TrackingManagerHelper::TrackNeutralParticle(aTrack, fStep, physics);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AdeptCPUTrackingManager::HandOverOneTrack(G4Track* aTrack)
{
    G4ParticleDefinition const* part = aTrack->GetParticleDefinition();

    if (part == G4Electron::Definition() || part == G4Positron::Definition())
    {
        TrackElectron(aTrack);
    }
    else if (part == G4Gamma::Definition())
    {
        TrackGamma(aTrack);
    }

    aTrack->SetTrackStatus(fStopAndKill);
    delete aTrack;
}
