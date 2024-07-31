#pragma once

#include "GPUOffloadInterface.hh"

//! Available offload backends
enum class GPUOffloadBackend
{
  Adept,
  Celeritas,
  CeleritasCPU,
  None,
  Unknown
};

//! Options for type of offload backend and its parameters
struct GPUOffloadOptions
{
  //! Implementation to use for offload backend
  GPUOffloadBackend backend = GPUOffloadBackend::Unknown;

  //! True if known backend selected
  explicit operator bool() const { return backend != GPUOffloadBackend::Unknown; }
};

//! Setup offload backend from input options
void GPUOffloadSetup(const GPUOffloadOptions& opts);

//! Get reference to offloader for the current thread
GPUOffloadInterface& GPUOffload();
