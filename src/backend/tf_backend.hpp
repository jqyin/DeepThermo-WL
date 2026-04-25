#pragma once

#include "backend/inference_backend.hpp"

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/core/protobuf/meta_graph.pb.h>
#include <tensorflow/core/public/session.h>
#include <tensorflow/core/public/session_options.h>

class TFBackend : public InferenceBackend {
public:
    void load(const std::string& dir, int local_gpu, int rank,
              const deepthermo::AlloyState& alloy) override;
    void encode(const float* input, int D, int NE, float* z) override;
    void decode(const float* z, int D, int NE, float* output) override;
    const char* name() const override { return "tensorflow"; }

private:
    tensorflow::SavedModelBundle encoder_bundle_;
    tensorflow::SavedModelBundle decoder_bundle_;
    void load_one(const std::string& path, int local_gpu,
                  tensorflow::SavedModelBundle& bundle);
    static void run_session(tensorflow::SavedModelBundle& bundle,
                            const std::string& input_name,
                            const tensorflow::Tensor& in,
                            tensorflow::Tensor& out);
};
