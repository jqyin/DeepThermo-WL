#pragma once

#include "backend/inference_backend.hpp"

#include <torch/script.h>

class TorchBackend : public InferenceBackend {
public:
    TorchBackend();
    ~TorchBackend() override;
    void load(const std::string& dir, int local_gpu, int rank,
              const deepthermo::AlloyState& alloy) override;
    void encode(const float* input, int D, int NE, float* z) override;
    void decode(const float* z, int D, int NE, float* output) override;
    const char* name() const override { return "libtorch"; }

private:
    static torch::jit::script::Module load_module(const std::string& dir,
                                                  const std::string& base,
                                                  const std::string& alloy_tag);
    torch::Device device_{torch::kCPU};
    torch::jit::script::Module encoder_;
    torch::jit::script::Module decoder_;
};
