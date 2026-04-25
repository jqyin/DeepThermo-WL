/**
 * @file smartredis_backend.hpp
 * @brief SmartRedis (RedisAI) implementation of ::InferenceBackend.
 */

#pragma once

#include "backend/inference_backend.hpp"

#include <memory>
#include <string>

#include "client.h"

/**
 * @brief Backend that talks to a co-located RedisAI server per rank.
 *
 * Used when configuring with `-DDEEPTHERMO_BACKEND=redis`. The launcher
 * (e.g. examples/perlmutter/hea-wl.sb under the legacy SmartRedis
 * configuration) is responsible for starting one redis-server per
 * GPU on every node. Each rank connects to `127.0.0.1:6781+local_gpu`.
 *
 * Models are shipped as frozen `.pb` files at `<dir>/encoder_<alloy>.pb`
 * and `<dir>/decoder_<alloy>.pb`.
 */
class SmartRedisBackend : public InferenceBackend {
public:
    /// @copydoc InferenceBackend::load
    void load(const std::string& dir, int local_gpu, int rank,
              const deepthermo::AlloyState& alloy) override;
    /// @copydoc InferenceBackend::encode
    void encode(const float* input, int D, int NE, float* z) override;
    /// @copydoc InferenceBackend::decode
    void decode(const float* z, int D, int NE, float* output) override;
    /// @copydoc InferenceBackend::name
    const char* name() const override { return "smartredis"; }

private:
    std::unique_ptr<SmartRedis::Client> client_; ///< Per-rank Redis client.
    std::string encoder_name_;                   ///< Server-side key prefix for encoder.
    std::string decoder_name_;                   ///< Server-side key prefix for decoder.
    std::string rank_str_;                       ///< Cached `to_string(rank)` for keys.
};
