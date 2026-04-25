/**
 * @file torch_backend.hpp
 * @brief LibTorch (TorchScript) implementation of ::InferenceBackend.
 */

#pragma once

#include "backend/inference_backend.hpp"

#include <torch/script.h>

/**
 * @brief Default (LibTorch) inference backend.
 *
 * Loads two TorchScript modules — `encoder_<alloy>.pt` and
 * `decoder_<alloy>.pt` — at simulation start and binds each MPI rank to
 * a local GPU. The same code path runs on CUDA (Perlmutter) and ROCm
 * (Frontier) since LibTorch exposes both under `torch::kCUDA`.
 *
 * The exported wrappers internally permute NHWDC↔NCDHW so the C++ side
 * never has to think about layout — see
 * `vae-modeling/vae/src/export_torchscript.py`.
 */
class TorchBackend : public InferenceBackend {
public:
    TorchBackend();
    ~TorchBackend() override;

    /// @copydoc InferenceBackend::load
    void load(const std::string& dir, int local_gpu, int rank,
              const deepthermo::AlloyState& alloy) override;
    /// @copydoc InferenceBackend::encode
    void encode(const float* input, int D, int NE, float* z) override;
    /// @copydoc InferenceBackend::decode
    void decode(const float* z, int D, int NE, float* output) override;
    /// @copydoc InferenceBackend::name
    const char* name() const override { return "libtorch"; }

private:
    /**
     * @brief Resolve `<dir>/<base>_<alloy_tag>.pt` (with bare `<base>.pt`
     *        as fallback) and load it.
     *
     * @throws std::runtime_error if no candidate exists.
     */
    static torch::jit::script::Module load_module(const std::string& dir,
                                                  const std::string& base,
                                                  const std::string& alloy_tag);

    torch::Device device_{torch::kCPU};      ///< Selected device.
    torch::jit::script::Module encoder_;     ///< Encoder TorchScript module.
    torch::jit::script::Module decoder_;     ///< Decoder TorchScript module.
};
