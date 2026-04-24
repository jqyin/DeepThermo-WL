// model.hpp — TensorFlow C++ SavedModel wrapper (legacy backend).
//
// Only included when DEEPTHERMO_BACKEND=tf. Task 4 of the v2 refactor will
// move this behind the InferenceBackend interface and into a proper
// translation unit; for now it's a header-only wrapper with `inline`
// applied so it can be safely included from multiple TUs without needing
// --allow-multiple-definition at link time.

#pragma once

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/core/protobuf/meta_graph.pb.h>
#include <tensorflow/core/public/session.h>
#include <tensorflow/core/public/session_options.h>

using tensorflow::RunOptions;
using tensorflow::SavedModelBundle;
using tensorflow::SessionOptions;
using tensorflow::Tensor;

enum PredictMode { encoder, decoder };

class heaModel {
public:
    void LoadModel(std::string model_path, SessionOptions options);
    void Predict(Tensor input, std::vector<Tensor>& pred, PredictMode mode);

private:
    SavedModelBundle bundle_;
    SessionOptions sessOpt_;
    RunOptions runOpt_;
};

inline void heaModel::LoadModel(std::string model_path, SessionOptions options) {
    sessOpt_ = options;
    sessOpt_.config.mutable_gpu_options()->set_allow_growth(true);
    auto status = tensorflow::LoadSavedModel(sessOpt_, runOpt_, model_path, {"serve"}, &bundle_);
    if (!status.ok()) {
        std::cout << "Error in loading model" << std::endl;
    }
}

inline void heaModel::Predict(Tensor input, std::vector<Tensor>& pred, PredictMode mode) {
    std::string input_node = (mode == encoder) ? "input_1" : "input_2";
    std::string output_node = "Identity";
    std::vector<std::pair<std::string, Tensor>> data = {{input_node, input}};
    auto status = bundle_.GetSession()->Run(data, {output_node}, {}, &pred);
    if (!status.ok()) {
        std::cout << "Error in predict" << std::endl;
        std::cout << status.error_message();
    }
}

// Shared encoder/decoder wrappers. Defined exactly once in main.cc via
// DEEPTHERMO_DEFINE_TF_MODELS.
inline std::array<heaModel, 2> tf_models;
