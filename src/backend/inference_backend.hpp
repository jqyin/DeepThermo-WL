/**
 * @file inference_backend.hpp
 * @brief Runtime-polymorphic encoder/decoder driver.
 *
 * Three concrete implementations exist, selected at build time by
 * `DEEPTHERMO_BACKEND` in CMake:
 *   - ::TorchBackend       — LibTorch (default)
 *   - ::TFBackend          — TensorFlow C++ SavedModel
 *   - ::SmartRedisBackend  — RedisAI via SmartRedis client
 *
 * Exactly one is compiled into `deepthermo_core`. ::make_backend
 * returns an instance of the selected one; callers don't need to know
 * which.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include "sim_context.hpp"

/**
 * @brief Pure-virtual interface for the VAE encoder/decoder pair.
 *
 * Implementations wrap whatever framework actually evaluates the model
 * (LibTorch / TensorFlow / SmartRedis) and adapt to/from the channels-
 * last (NHWDC) buffer layout that the simulation core produces.
 */
class InferenceBackend {
public:
    virtual ~InferenceBackend() = default;

    /**
     * @brief Set up the backend and load the encoder + decoder.
     *
     * @param dir       Directory containing the model files (TF SavedModel
     *                  sub-dirs, frozen `.pb`, or TorchScript `.pt`).
     * @param local_gpu Which device on this node to bind to.
     * @param rank      Global MPI rank (used by per-rank-keyed backends
     *                  such as SmartRedis).
     * @param alloy     Read-only access to `NE` and the element table;
     *                  used for filename composition (e.g.
     *                  `encoder_MoNbTaW.pt`).
     *
     * @throws std::runtime_error on missing model files or framework
     *         load failure.
     */
    virtual void load(const std::string& dir, int local_gpu, int rank,
                      const deepthermo::AlloyState& alloy) = 0;

    /**
     * @brief Run the encoder forward pass.
     *
     * @param[in]  input D*D*D*NE one-hot floats (NHWDC layout).
     * @param[in]  D     Spatial grid dimension (= `alloyState.VAE_D`).
     * @param[in]  NE    Number of element channels.
     * @param[out] z     Length-3 latent vector.
     */
    virtual void encode(const float* input, int D, int NE, float* z) = 0;

    /**
     * @brief Run the decoder forward pass.
     *
     * @param[in]  z      Length-3 latent vector.
     * @param[in]  D      Spatial grid dimension.
     * @param[in]  NE     Number of element channels.
     * @param[out] output D*D*D*NE per-site, per-element probabilities (NHWDC).
     */
    virtual void decode(const float* z, int D, int NE, float* output) = 0;

    /**
     * @brief Identifier used in log lines and error messages.
     */
    virtual const char* name() const = 0;
};

/**
 * @brief Return an instance of the build-selected InferenceBackend.
 *
 * Defined in `backend/backend_factory.cc`; resolves at compile time
 * based on the `DEEPTHERMO_BACKEND_TORCH` / `_TF` / `_REDIS` macro.
 */
std::unique_ptr<InferenceBackend> make_backend();
