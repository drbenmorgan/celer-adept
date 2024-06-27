# Offloading to GPU using G4VTrackingManager

Testbed application for offloading from Geant4 to AdePT/Celeritas via
concrete `G4VTrackingManager`s.

- Understand commonalities/differences between how problem data is setup
- What setup/offload options are common
- Future common testing problems (easy switching between impls)

It's essentially (intended to be) a merge/blend of Celeritas' `accel` 
examples and AdePT's `example1`. At present only a trivial setup is
used to focus on the off/onload.

# Build/Use
The project uses CMake's `FetchContent` to pull in the relevant branches
of AdePT, Celeritas, and the supporting G4VG Geant4-VecGeom in memory converter.
WIth these, the build should then work with either AdePT's LCG view or a 
Celeritas Spack view (if not, then relevant patches here or to the upstreams
should be submitted).

**NB**: One current issue is that the _first_ invocation of CMake will result in
an error:

```console
-- Performing Test Geant4_flush_FOUND
CMake Error in /tmp/adept-build/ca-build/CMakeFiles/CMakeScratch/TryCompile-7OCbcN/CMakeLists.txt:
  No known features for CUDA compiler

  ""

  version .


CMake Error at /usr/share/cmake/Modules/Internal/CheckSourceCompiles.cmake:101 (try_compile):
  Failed to generate test project build system.
Call Stack (most recent call first):
  /usr/share/cmake/Modules/CheckCXXSourceCompiles.cmake:76 (cmake_check_source_compiles)
  /tmp/adept-build/ca-build/_deps/adept-src/CMakeLists.txt:105 (check_cxx_source_compiles)
```

**Just running `cmake` again will resolve the issue and yield a successful configuration,
and then build.**

The only program is the `trackingmanager-offload`, so just run this to
test at the moment. Comments in `trackingmanager-offload.cc` should illustrate
what's going on (and provide a basis for Doxygen going forward).

This hardcodes the use of Celeritas offload, but this can be changed for
the AdePT and no offload cases in `GPUOffload/GPUOffload.hh`.

# Known Limitations
## Only one Run is permitted after initialization
As shown in [`trackingmanager-offload.cc`](./trackingmanager-offload.cc), only one
call to `G4RunManager::BeamOn` can be made per application execution. This is due
to the potential for geometry/physics to change between runs which would require
reinitialization of GPU-side data for AdePT/Celeritas. Neither project supports
this capability yet, as the primary HEP production use case does not require this.
This may be added later, but for now only one run is supported.
