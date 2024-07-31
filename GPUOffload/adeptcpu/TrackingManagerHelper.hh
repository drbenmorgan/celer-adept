//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// TrackingManagerHelper (copied from examples/extended/runAndEvent/RE07)
//
// Class description:
//
// Helper class for reducing the effort required to implement a custom tracking
// manager. It implements a stepping loop that calls user actions as the
// generic tracking and stepping managers do, and it implements navigation for
// charged particles in energy-preserving fields and for neutral particles.
//
// Original author: Jonas Hahnfeld, 2021

#ifndef TrackingManagerHelper_hh
#define TrackingManagerHelper_hh 1

#include <G4EventManager.hh>
#include <G4Field.hh>
#include <G4FieldManager.hh>
#include <G4FieldManagerStore.hh>
#include <G4GeometryTolerance.hh>
#include <G4LogicalVolume.hh>
#include <G4Navigator.hh>
#include <G4PropagatorInField.hh>
#include <G4Region.hh>
#include <G4SafetyHelper.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4ThreeVector.hh>
#include <G4TouchableHandle.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4TrackVector.hh>
#include <G4TransportationManager.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSensitiveDetector.hh>
#include <globals.hh>

class TrackingManagerHelper
{
  public:
    class Physics
    {
      public:
        virtual void StartTracking(G4Track*) {}
        virtual void EndTracking() {}

        // Combines AlongStep and PostStep; the implementation needs to
        // remember the right value to pass as previousStepSize to G4VProcess.
        virtual G4double GetPhysicalInteractionLength(G4Track const& track) = 0;

        // This method is called for every step after navigation. The updated
        // position is stored in the G4Step's post-step point. Any particle
        // change should be applied directly to the step, UpdateTrack() will be
        // called automatically after this method returns. If secondaries
        // should be given back to the G4EventManager, put them into the
        // container passed as the last argument.
        virtual void
        AlongStepDoIt(G4Track& track, G4Step& step, G4TrackVector& secondaries)
            = 0;

        // This method is called unless the track has been killed during this
        // step. If secondaries should be given back to the G4EventManager, put
        // them into the container passed as the last argument.
        virtual void
        PostStepDoIt(G4Track& track, G4Step& step, G4TrackVector& secondaries)
            = 0;

        virtual bool HasAtRestProcesses() { return false; }

        // This method is called when a track is stopped, but still alive. If
        // secondaries should be given back to the G4EventManager, put them
        // into the container passed as the last argument.
        virtual void
        AtRestDoIt(G4Track& track, G4Step& step, G4TrackVector& secondaries)
        {
            (void)track;
            (void)step;
            (void)secondaries;
        }
    };

    class Navigation
    {
      public:
        virtual G4double
        MakeStep(G4Track& track, G4Step& step, G4double physicalStep)
            = 0;

        virtual void FinishStep(G4Track& track, G4Step& step) = 0;
    };

    class ChargedNavigation final : public Navigation
    {
      public:
        inline ChargedNavigation();
        inline G4double
        MakeStep(G4Track& track, G4Step& step, G4double physicalStep) override;
        inline void FinishStep(G4Track& track, G4Step& step) override;

      private:
        G4Navigator* fLinearNavigator;
        G4PropagatorInField* fFieldPropagator;
        G4SafetyHelper* fSafetyHelper;
        G4ThreeVector fSafetyOrigin;
        G4double fSafety = 0;
        G4double fPostStepSafety = 0;
        G4double kCarTolerance;
        G4bool fGeometryLimitedStep;
    };

    class NeutralNavigation final : public Navigation
    {
      public:
        inline NeutralNavigation();
        inline G4double
        MakeStep(G4Track& track, G4Step& step, G4double physicalStep) override;
        inline void FinishStep(G4Track& track, G4Step& step) override;

      private:
        G4Navigator* fLinearNavigator;
        G4SafetyHelper* fSafetyHelper;
        G4ThreeVector fSafetyOrigin;
        G4double fSafety = 0;
        G4double fPostStepSafety = 0;
        G4double kCarTolerance;
        G4bool fGeometryLimitedStep;
    };

    template<typename PhysicsImpl, typename NavigationImpl>
    static void TrackParticle(G4Track* aTrack,
                              G4Step* aStep,
                              PhysicsImpl& physics,
                              NavigationImpl& navigation);

    template<typename PhysicsImpl>
    static void
    TrackChargedParticle(G4Track* aTrack, G4Step* aStep, PhysicsImpl& physics);

    template<typename PhysicsImpl>
    static void
    TrackNeutralParticle(G4Track* aTrack, G4Step* aStep, PhysicsImpl& physics);
};

template<typename PhysicsImpl, typename NavigationImpl>
void TrackingManagerHelper::TrackParticle(G4Track* aTrack,
                                          G4Step* aStep,
                                          PhysicsImpl& physics,
                                          NavigationImpl& navigation)
{
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
    G4Step& step = *aStep;
    G4TrackVector& secondaries = *step.GetfSecondary();
    G4StepPoint& preStepPoint = *step.GetPreStepPoint();
    step.InitializeStep(aTrack);
    aTrack->SetStep(&step);

    // Start of tracking: Inform user and processes.
    if (userTrackingAction)
    {
        userTrackingAction->PreUserTrackingAction(aTrack);
    }

    physics.StartTracking(aTrack);

    while (aTrack->GetTrackStatus() == fAlive)
    {
        // Beginning of this step: Prepare data structures.
        aTrack->IncrementCurrentStepNumber();

        step.CopyPostToPreStepPoint();
        step.ResetTotalEnergyDeposit();
        aTrack->SetTouchableHandle(aTrack->GetNextTouchableHandle());

        auto* lvol = aTrack->GetTouchable()->GetVolume()->GetLogicalVolume();
        preStepPoint.SetMaterial(lvol->GetMaterial());
        preStepPoint.SetMaterialCutsCouple(lvol->GetMaterialCutsCouple());

        // Query step lengths from pyhsics and geometry, decide on limit.
        G4double physicalStep = physics.GetPhysicalInteractionLength(*aTrack);
        G4double geometryStep
            = navigation.MakeStep(*aTrack, step, physicalStep);

        bool geometryLimitedStep = geometryStep < physicalStep;
        G4double finalStep = geometryLimitedStep ? geometryStep : physicalStep;

        step.SetStepLength(finalStep);
        aTrack->SetStepLength(finalStep);

        // Call AlongStepDoIt in every step.
        physics.AlongStepDoIt(*aTrack, step, secondaries);
        step.UpdateTrack();

        if (aTrack->GetTrackStatus() == fAlive
            && aTrack->GetKineticEnergy() < DBL_MIN)
        {
            if (physics.HasAtRestProcesses())
            {
                aTrack->SetTrackStatus(fStopButAlive);
            }
            else
            {
                aTrack->SetTrackStatus(fStopAndKill);
            }
        }

        navigation.FinishStep(*aTrack, step);

        // Check if the track left the world.
        if (aTrack->GetNextVolume() == nullptr)
        {
            aTrack->SetTrackStatus(fStopAndKill);
        }

        // The check should rather check for == fAlive and avoid calling
        // PostStepDoIt for fStopButAlive, but the generic stepping loop
        // does it like this...
        if (aTrack->GetTrackStatus() != fStopAndKill)
        {
            physics.PostStepDoIt(*aTrack, step, secondaries);
        }

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

    if (aTrack->GetTrackStatus() == fStopButAlive
        && aTrack->GetNextVolume() != nullptr)
    {
        // Do one final step.
        aTrack->IncrementCurrentStepNumber();

        step.CopyPostToPreStepPoint();
        step.ResetTotalEnergyDeposit();

        physics.AtRestDoIt(*aTrack, step, secondaries);

        // End of this step: Call sensitive detector and stepping actions.
        auto* lvol = aTrack->GetTouchable()->GetVolume()->GetLogicalVolume();
        if (step.GetControlFlag() != AvoidHitInvocation)
        {
            auto sensitive = lvol->GetSensitiveDetector();
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
    physics.EndTracking();

    if (userTrackingAction)
    {
        userTrackingAction->PostUserTrackingAction(aTrack);
    }

    evtMgr->StackTracks(&secondaries);
}

TrackingManagerHelper::ChargedNavigation::ChargedNavigation()
{
    auto* transMgr = G4TransportationManager::GetTransportationManager();
    fLinearNavigator = transMgr->GetNavigatorForTracking();
    fFieldPropagator = transMgr->GetPropagatorInField();
    fSafetyHelper = transMgr->GetSafetyHelper();
    kCarTolerance
        = 0.5 * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

    // Reset sstate of field propagator and all chord finders.
    fFieldPropagator->ClearPropagatorState();

    auto* fieldMgrStore = G4FieldManagerStore::GetInstance();
    fieldMgrStore->ClearAllChordFindersState();
}

G4double
TrackingManagerHelper::ChargedNavigation::MakeStep(G4Track& track,
                                                   G4Step& step,
                                                   G4double physicalStep)
{
    G4ThreeVector pos = track.GetPosition();
    G4ThreeVector dir = track.GetMomentumDirection();
    G4StepPoint& postStepPoint = *step.GetPostStepPoint();

    bool fieldExertsForce = false;
    if (auto* fieldMgr
        = fFieldPropagator->FindAndSetFieldManager(track.GetVolume()))
    {
        fieldMgr->ConfigureForTrack(&track);
        if (G4Field const* ptrField = fieldMgr->GetDetectorField())
        {
            fieldExertsForce = true;
        }
    }

    G4double endpointDistance;
    G4double safety = 0.0;
    // Setting a fallback value for safety is required in case of where very
    // short steps where the field propagator returns immediately without
    // calling geometry.
    G4double const shiftSquare = (pos - fSafetyOrigin).mag2();
    if (shiftSquare < sqr(fSafety))
    {
        safety = fSafety - std::sqrt(shiftSquare);
    }

    if (fieldExertsForce)
    {
        G4DynamicParticle const* pParticle = track.GetDynamicParticle();
        G4double const particleCharge = pParticle->GetCharge();
        G4double const particleMass = pParticle->GetMass();
        G4double const magneticMoment = pParticle->GetMagneticMoment();
        G4ThreeVector const particleSpin = pParticle->GetPolarization();
        G4double const kineticEnergy = pParticle->GetKineticEnergy();
        auto const pParticleDef = pParticle->GetDefinition();
        auto const particlePDGSpin = pParticleDef->GetPDGSpin();
        auto const particlePDGMagM = pParticleDef->GetPDGMagneticMoment();

        auto equationOfMotion = fFieldPropagator->GetCurrentEquationOfMotion();
        equationOfMotion->SetChargeMomentumMass(
            G4ChargeState(particleCharge, magneticMoment, particlePDGSpin),
            pParticle->GetTotalMomentum(),
            particleMass);

        G4ThreeVector const startPosition = pos;
        G4ThreeVector const startDirection = dir;
        G4FieldTrack aFieldTrack(startPosition,
                                 track.GetGlobalTime(),  // Lab.
                                 dir,
                                 kineticEnergy,
                                 particleMass,
                                 particleCharge,
                                 particleSpin,
                                 particlePDGMagM,
                                 0.0,  // Length along track
                                 particlePDGSpin);

        // Do the Transport in the field (non recti-linear)
        //
        fGeometryLimitedStep = false;
        G4double const lengthAlongCurve
            = fFieldPropagator->ComputeStep(aFieldTrack,
                                            physicalStep,
                                            safety,
                                            track.GetVolume(),
                                            kineticEnergy < 250.0);
        if (lengthAlongCurve < physicalStep)
        {
            physicalStep = lengthAlongCurve;
            fGeometryLimitedStep = true;
        }
        fSafetyHelper->SetCurrentSafety(safety, pos);
        fSafetyOrigin = pos;
        fSafety = safety;

        if (fFieldPropagator->IsParticleLooping())
        {
            track.SetTrackStatus(fStopAndKill);
        }

        pos = aFieldTrack.GetPosition();
        dir = aFieldTrack.GetMomentumDir();

        postStepPoint.SetPosition(pos);
        postStepPoint.SetMomentumDirection(dir);

        endpointDistance = (startPosition - pos).mag();
    }
    else
    {
        fGeometryLimitedStep = false;
        G4double linearStepLength
            = fLinearNavigator->ComputeStep(pos, dir, physicalStep, safety);
        if (linearStepLength < physicalStep)
        {
            physicalStep = linearStepLength;
            fGeometryLimitedStep = true;
        }
        fSafetyHelper->SetCurrentSafety(safety, pos);
        fSafetyOrigin = pos;
        fSafety = safety;

        // Update the position.
        pos += physicalStep * dir;
        postStepPoint.SetPosition(pos);

        endpointDistance = physicalStep;
    }

    // Update global, local, and proper time.
    double velocity = track.GetVelocity();
    double deltaTime = 0;
    if (velocity > 0)
    {
        deltaTime = physicalStep / velocity;
    }

    postStepPoint.AddGlobalTime(deltaTime);
    postStepPoint.AddLocalTime(deltaTime);

    double restMass = track.GetDynamicParticle()->GetMass();
    double deltaProperTime = deltaTime * (restMass / track.GetTotalEnergy());
    postStepPoint.AddProperTime(deltaProperTime);

    // Compute safety, including the call to safetyHelper, but don't set the
    // safety in the post-step point to mimick the generic stepping loop.
    if (safety > physicalStep)
    {
        safety -= physicalStep;
    }
    else if (safety < endpointDistance)
    {
        safety = fLinearNavigator->ComputeSafety(pos);
        fSafetyHelper->SetCurrentSafety(safety, pos);
        fSafetyOrigin = pos;
        fSafety = safety;
    }
    else
    {
        safety = 0;
    }
    if (safety < kCarTolerance)
    {
        fPostStepSafety = kCarTolerance;
    }
    else
    {
        fPostStepSafety = safety;
    }

    return physicalStep;
}

void TrackingManagerHelper::ChargedNavigation::FinishStep(G4Track& track,
                                                          G4Step& step)
{
    // Now set the safety that was computed in MakeStep.
    G4StepPoint& postStepPoint = *step.GetPostStepPoint();
    postStepPoint.SetSafety(fPostStepSafety);

    G4TouchableHandle touchableHandle = track.GetTouchableHandle();
    G4ThreeVector const& pos = track.GetPosition();
    if (fGeometryLimitedStep)
    {
        // Relocate the particle.
        fLinearNavigator->SetGeometricallyLimitedStep();
        fLinearNavigator->LocateGlobalPointAndUpdateTouchableHandle(
            pos, track.GetMomentumDirection(), touchableHandle, true);
        G4VPhysicalVolume const* newVolume = touchableHandle->GetVolume();
        if (newVolume == nullptr)
        {
            postStepPoint.SetStepStatus(fWorldBoundary);
        }
        else
        {
            postStepPoint.SetStepStatus(fGeomBoundary);
        }
    }
    else
    {
        // Move the Navigator's location.
        fLinearNavigator->LocateGlobalPointWithinVolume(pos);
    }

    postStepPoint.SetTouchableHandle(touchableHandle);
    track.SetNextTouchableHandle(touchableHandle);
}

template<typename PhysicsImpl>
void TrackingManagerHelper::TrackChargedParticle(G4Track* aTrack,
                                                 G4Step* aStep,
                                                 PhysicsImpl& physics)
{
    ChargedNavigation navigation;
    TrackParticle(aTrack, aStep, physics, navigation);
}

TrackingManagerHelper::NeutralNavigation::NeutralNavigation()
{
    auto* transMgr = G4TransportationManager::GetTransportationManager();
    fLinearNavigator = transMgr->GetNavigatorForTracking();
    fSafetyHelper = transMgr->GetSafetyHelper();
    kCarTolerance
        = 0.5 * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

G4double
TrackingManagerHelper::NeutralNavigation::MakeStep(G4Track& track,
                                                   G4Step& step,
                                                   G4double physicalStep)
{
    G4ThreeVector pos = track.GetPosition();
    G4ThreeVector dir = track.GetMomentumDirection();
    G4StepPoint& postStepPoint = *step.GetPostStepPoint();

    G4double safety = 0.0;
    G4double const shiftSquare = (pos - fSafetyOrigin).mag2();
    if (shiftSquare < sqr(fSafety))
    {
        safety = fSafety - std::sqrt(shiftSquare);
    }

    fGeometryLimitedStep = false;
    G4double linearStepLength
        = fLinearNavigator->ComputeStep(pos, dir, physicalStep, safety);
    if (linearStepLength < physicalStep)
    {
        physicalStep = linearStepLength;
        fGeometryLimitedStep = true;
    }
    fSafetyHelper->SetCurrentSafety(safety, pos);
    fSafetyOrigin = pos;
    fSafety = safety;

    // Update the position.
    pos += physicalStep * dir;
    postStepPoint.SetPosition(pos);

    // Update global, local, and proper time.
    double velocity = track.GetVelocity();
    double deltaTime = 0;
    if (velocity > 0)
    {
        deltaTime = physicalStep / velocity;
    }
    postStepPoint.AddGlobalTime(deltaTime);
    postStepPoint.AddLocalTime(deltaTime);

    double restMass = track.GetDynamicParticle()->GetMass();
    double deltaProperTime = deltaTime * (restMass / track.GetTotalEnergy());
    postStepPoint.AddProperTime(deltaProperTime);

    // Compute safety, but don't set the safety in the post-step point to
    // mimick the generic stepping loop.
    if (safety > physicalStep)
    {
        safety -= physicalStep;
    }
    else
    {
        safety = 0;
    }
    if (safety < kCarTolerance)
    {
        fPostStepSafety = kCarTolerance;
    }
    else
    {
        fPostStepSafety = safety;
    }

    return physicalStep;
}

void TrackingManagerHelper::NeutralNavigation::FinishStep(G4Track& track,
                                                          G4Step& step)
{
    // Now set the safety that was computed in MakeStep.
    G4StepPoint& postStepPoint = *step.GetPostStepPoint();
    postStepPoint.SetSafety(fPostStepSafety);

    G4TouchableHandle touchableHandle = track.GetTouchableHandle();
    G4ThreeVector const& pos = track.GetPosition();
    if (fGeometryLimitedStep)
    {
        // Relocate the particle.
        fLinearNavigator->SetGeometricallyLimitedStep();
        fLinearNavigator->LocateGlobalPointAndUpdateTouchableHandle(
            pos, track.GetMomentumDirection(), touchableHandle, true);
        G4VPhysicalVolume const* newVolume = touchableHandle->GetVolume();
        if (newVolume == nullptr)
        {
            postStepPoint.SetStepStatus(fWorldBoundary);
        }
        else
        {
            postStepPoint.SetStepStatus(fGeomBoundary);
        }
    }
    else
    {
        // Move the Navigator's location.
        fLinearNavigator->LocateGlobalPointWithinVolume(pos);
    }

    postStepPoint.SetTouchableHandle(touchableHandle);
    track.SetNextTouchableHandle(touchableHandle);
}

template<typename PhysicsImpl>
void TrackingManagerHelper::TrackNeutralParticle(G4Track* aTrack,
                                                 G4Step* aStep,
                                                 PhysicsImpl& physics)
{
    NeutralNavigation navigation;
    TrackParticle(aTrack, aStep, physics, navigation);
}

#endif
