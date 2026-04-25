/**
 * @file tf_backend.hpp
 * @brief TensorFlow C++ SavedModel implementation of ::InferenceBackend.
 */

#pragma once

#include "backend/inference_backend.hpp"

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/core/protobuf/meta_graph.pb.h>
#include <tensorflow/core/public/session.h>
#include <tensorflow/core/public/session_options.h>

/**
 * @brief Backend that loads two TF SavedModels at `<dir>/encoder` and
 *        `<dir>/decoder`.
 *
 * Used when configuring with `-DDEEPTHERMO_BACKEND=tf`. Expected graph
 * I/O nodes are `input_1` (encoder) / `input_2` (decoder) and `Identity`
 * (output). Retraining must preserve those names.
 */
class TFBackend : public InferenceBackend {
public:
    /// @copydoc InferenceBackend::load
    void load(const std::string& dir, int local_gpu, int rank,
              const deepthermo::AlloyState& alloy) override;
    /// @copydoc InferenceBackend::encode
    void encode(const float* input, int D, int NE, float* z) override;
    /// @copydoc InferenceBackend::decode
    void decode(const float* z, int D, int NE, float* output) override;
    /// @copydoc InferenceBackend::name
    const char* name() const override { return "tensorflow"; }

private:
    tensorflow::SavedModelBundle encoder_bundle_; ///< Loaded encoder graph.
    tensorflow::SavedModelBundle decoder_bundle_; ///< Loaded decoder graph.

    /// @brief Helper used by ::load to load one SavedModel onto `local_gpu`.
    void load_one(const std::string& path, int local_gpu,
                  tensorflow::SavedModelBundle& bundle);
    /// @brief Helper used by ::encode / ::decode to invoke a session.
    static void run_session(tensorflow::SavedModelBundle& bundle,
                            const std::string& input_name,
                            const tensorflow::Tensor& in,
                            tensorflow::Tensor& out);
};
