#ifndef AdeptCPUTrackingManager_h
#define AdeptCPUTrackingManager_h 1

#include "G4VTrackingManager.hh"
#include "globals.hh"

class G4HepEmRunManager;
class G4HepEmRandomEngine;
class G4HepEmNoProcess;
class G4SafetyHelper;
class G4Step;

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AdeptCPUTrackingManager : public G4VTrackingManager {
public:
  AdeptCPUTrackingManager();
  ~AdeptCPUTrackingManager();

  void BuildPhysicsTable(const G4ParticleDefinition &) override;

  void PreparePhysicsTable(const G4ParticleDefinition &) override;

  void HandOverOneTrack(G4Track *aTrack) override;

  void SetMultipleSteps(G4bool val) {
    fMultipleSteps = val;
  }
  G4bool MultipleSteps() const {
    return fMultipleSteps;
  }

private:
  void TrackElectron(G4Track *aTrack);
  void TrackGamma(G4Track *aTrack);

  G4HepEmRunManager *fRunManager;
  G4HepEmRandomEngine *fRandomEngine;
  G4SafetyHelper *fSafetyHelper;
  G4Step *fStep;

  const std::vector<G4double> *theCutsGamma = nullptr;
  const std::vector<G4double> *theCutsElectron = nullptr;
  const std::vector<G4double> *theCutsPositron = nullptr;
  G4bool applyCuts = false;
  G4bool fMultipleSteps = true;

  // A set of empty processes with the correct names and types just to be able
  // to set them as process limiting the step and creating secondaries as some
  // user codes rely on this information.
  std::vector<G4HepEmNoProcess *> fElectronNoProcessVector;
  std::vector<G4HepEmNoProcess *> fGammaNoProcessVector;
  G4HepEmNoProcess *fTransportNoProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
