#include "backend/tf_backend.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

void TFBackend::load_one(const std::string& path, int local_gpu,
                         tensorflow::SavedModelBundle& bundle) {
    tensorflow::SessionOptions opt;
    opt.config.mutable_gpu_options()->set_visible_device_list(std::to_string(local_gpu));
    opt.config.mutable_gpu_options()->set_allow_growth(true);
    tensorflow::RunOptions run_opt;
    auto status = tensorflow::LoadSavedModel(opt, run_opt, path, {"serve"}, &bundle);
    if (!status.ok()) {
        throw std::runtime_error("TFBackend: failed to load " + path);
    }
}

void TFBackend::load(const std::string& dir, int local_gpu, int /*rank*/,
                     const deepthermo::AlloyState& /*alloy*/) {
    load_one(dir + "/encoder", local_gpu, encoder_bundle_);
    load_one(dir + "/decoder", local_gpu, decoder_bundle_);
}

void TFBackend::run_session(tensorflow::SavedModelBundle& bundle,
                            const std::string& input_name,
                            const tensorflow::Tensor& in,
                            tensorflow::Tensor& out) {
    std::vector<std::pair<std::string, tensorflow::Tensor>> feed = {{input_name, in}};
    std::vector<tensorflow::Tensor> outs;
    auto status = bundle.GetSession()->Run(feed, {"Identity"}, {}, &outs);
    if (!status.ok()) {
        throw std::runtime_error(std::string{"TFBackend: session run failed: "} +
                                 status.error_message());
    }
    out = std::move(outs[0]);
}

void TFBackend::encode(const float* input, int D, int NE, float* z) {
    tensorflow::TensorShape shape({1, D, D, D, NE});
    tensorflow::Tensor in(tensorflow::DT_FLOAT, shape);
    std::memcpy(in.flat<float>().data(), input, sizeof(float) * D * D * D * NE);
    tensorflow::Tensor out;
    run_session(encoder_bundle_, "input_1", in, out);
    std::memcpy(z, out.flat<float>().data(), sizeof(float) * 3);
}

void TFBackend::decode(const float* z, int D, int NE, float* output) {
    tensorflow::TensorShape shape({1, 3});
    tensorflow::Tensor in(tensorflow::DT_FLOAT, shape);
    std::memcpy(in.flat<float>().data(), z, sizeof(float) * 3);
    tensorflow::Tensor out;
    run_session(decoder_bundle_, "input_2", in, out);
    std::memcpy(output, out.flat<float>().data(),
                sizeof(float) * static_cast<std::size_t>(D) * D * D * NE);
}
