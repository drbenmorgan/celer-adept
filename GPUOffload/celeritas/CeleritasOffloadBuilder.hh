#pragma once

#include <memory>

#include "GPUOffload.hh"

#include "CeleritasOffload.hh"

std::unique_ptr<CeleritasOffload>
BuildCeleritasOffload(GPUOffloadOptions const& opts);
