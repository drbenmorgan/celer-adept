# Offloading to GPU using G4VTrackingManager

Testbed application for offloading from Geant4 to AdePT/Celeritas via
concrete `G4VTrackingManager`s.

- Understand commonalities/differences between how problem data is setup
- What setup/offload options are common
- Future common testing problems (easy switching between impls)

It's essentially (intended to be) a merge/blend of Celeritas' `accel` 
examples and AdePT's `Example22`. At present only a trivial setup is
used.

# Build/Use
The project uses CMake's `FetchContent` to pull in AdePT's `master`
branch and Celeritas' v0.4.1 tag. The build should then work with
either AdePT's LCG view or a Celeritas Spack view (if not, then relevant
patches here or to the upstreams should be submitted).

The only program is the `trackingmanager-offload`, so just run this to
test at the moment. Comments in `trackingmanager-offload.cc` should illustrate
what's going on (and provide a basis for Doxygen going forward).
