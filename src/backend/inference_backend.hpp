// inference_backend.hpp — runtime-polymorphic encoder/decoder driver.
//
// Three concrete implementations exist:
//   torch_backend       — LibTorch (default)
//   tf_backend          — TensorFlow C++ SavedModel
//   smartredis_backend  — RedisAI via SmartRedis client
//
// Exactly one is selected at build time by CMake (DEEPTHERMO_BACKEND=...);
// make_backend() returns an instance of the chosen one. Callers don't need
// to know which.

#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include "sim_context.hpp"

class InferenceBackend {
public:
    virtual ~InferenceBackend() = default;

    // Set up the backend: load encoder + decoder models, bind a GPU, etc.
    //   dir       — directory containing the model files (TF SavedModel
    //               sub-dirs, frozen .pb, or TorchScript .pt files)
    //   local_gpu — which device on this node to bind to
    //   rank      — global MPI rank (used by per-rank-keyed backends)
    //   alloy     — read-only access to NE and the element table; used
    //               for naming (e.g. "encoder_MoNbTaW.pt")
    virtual void load(const std::string& dir, int local_gpu, int rank,
                      const deepthermo::AlloyState& alloy) = 0;

    // Encode the current configuration into the latent space.
    // input:  D*D*D*NE one-hot floats (NHWDC layout)
    // z:      length-3 latent vector
    virtual void encode(const float* input, int D, int NE, float* z) = 0;

    // Decode a latent vector back to per-site probabilities.
    // z:      length-3 latent vector
    // output: D*D*D*NE floats (per-site, per-element probability)
    virtual void decode(const float* z, int D, int NE, float* output) = 0;

    virtual const char* name() const = 0;
};

// Factory: returns the compile-time selected backend.
std::unique_ptr<InferenceBackend> make_backend();
