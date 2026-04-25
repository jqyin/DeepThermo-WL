// torch_backend.cc — LibTorch (TorchScript) inference backend.
//
// Loads two TorchScript modules (encoder and decoder) at simulation start
// and binds each MPI rank to one local GPU. The same code path runs on
// CUDA (Perlmutter A100) and ROCm (Frontier MI250X) — LibTorch exposes
// both under torch::kCUDA, so there's no per-platform compile branch.

#include "backend/torch_backend.hpp"

#include <torch/cuda.h>

#include <cstring>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

TorchBackend::TorchBackend() = default;
TorchBackend::~TorchBackend() = default;

torch::jit::script::Module TorchBackend::load_module(const std::string& dir,
                                                     const std::string& base,
                                                     const std::string& alloy_tag) {
    // Prefer the alloy-tagged file (e.g. encoder_MoNbTaW.pt); fall back to
    // plain encoder.pt for tests / single-alloy installs.
    const std::vector<fs::path> candidates = {
        fs::path(dir) / (base + "_" + alloy_tag + ".pt"),
        fs::path(dir) / (base + ".pt"),
    };
    for (const auto& path : candidates) {
        if (fs::exists(path)) {
            try {
                return torch::jit::load(path.string());
            } catch (const c10::Error& e) {
                throw std::runtime_error("TorchBackend: failed to load " + path.string() +
                                         ": " + e.what());
            }
        }
    }
    throw std::runtime_error("TorchBackend: no model file found for '" + base +
                             "' in " + dir);
}

void TorchBackend::load(const std::string& dir, int local_gpu, int rank,
                        const deepthermo::AlloyState& alloy) {
    if (torch::cuda::is_available()) {
        device_ = torch::Device(torch::kCUDA, local_gpu);
    } else {
        device_ = torch::Device(torch::kCPU);
        if (rank == 0) {
            std::cerr << "TorchBackend: CUDA unavailable, falling back to CPU."
                      << std::endl;
        }
    }

    std::string alloy_tag;
    for (int i = 0; i < alloy.NE; ++i) alloy_tag += alloy.element[i + 1];

    encoder_ = load_module(dir, "encoder", alloy_tag);
    encoder_.to(device_);
    encoder_.eval();

    decoder_ = load_module(dir, "decoder", alloy_tag);
    decoder_.to(device_);
    decoder_.eval();

    if (rank == 0) {
        std::cerr << "TorchBackend: encoder + decoder loaded on " << device_
                  << std::endl;
    }
}

void TorchBackend::encode(const float* input, int D, int NE, float* z) {
    torch::NoGradGuard nograd;
    auto opts = torch::TensorOptions().dtype(torch::kFloat32);
    auto in = torch::from_blob(const_cast<float*>(input),
                               {1, D, D, D, NE}, opts).clone().to(device_);
    auto out = encoder_.forward({in}).toTensor().to(torch::kCPU).contiguous();
    if (out.numel() < 3) {
        throw std::runtime_error("TorchBackend::encode: encoder output has < 3 elements");
    }
    std::memcpy(z, out.data_ptr<float>(), 3 * sizeof(float));
}

void TorchBackend::decode(const float* z, int D, int NE, float* output) {
    torch::NoGradGuard nograd;
    auto opts = torch::TensorOptions().dtype(torch::kFloat32);
    auto in = torch::from_blob(const_cast<float*>(z),
                               {1, 3}, opts).clone().to(device_);
    auto out = decoder_.forward({in}).toTensor().to(torch::kCPU).contiguous();
    const std::size_t n = static_cast<std::size_t>(D) * D * D * NE;
    if (static_cast<std::size_t>(out.numel()) < n) {
        throw std::runtime_error("TorchBackend::decode: decoder output too small");
    }
    std::memcpy(output, out.data_ptr<float>(), n * sizeof(float));
}
