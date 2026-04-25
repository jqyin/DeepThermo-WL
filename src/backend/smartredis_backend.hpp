#pragma once

#include "backend/inference_backend.hpp"

#include <memory>
#include <string>

#include "client.h"

class SmartRedisBackend : public InferenceBackend {
public:
    void load(const std::string& dir, int local_gpu, int rank,
              const deepthermo::AlloyState& alloy) override;
    void encode(const float* input, int D, int NE, float* z) override;
    void decode(const float* z, int D, int NE, float* output) override;
    const char* name() const override { return "smartredis"; }

private:
    std::unique_ptr<SmartRedis::Client> client_;
    std::string encoder_name_;
    std::string decoder_name_;
    std::string rank_str_;
};
