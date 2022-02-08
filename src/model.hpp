#ifndef  MODEL_H
#define MODEL_H
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <tensorflow/core/public/session.h>
#include <tensorflow/core/protobuf/meta_graph.pb.h>
#include <tensorflow/core/public/session_options.h>
#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/cc/saved_model/tag_constants.h>

using namespace std;
//using tensorflow::string;
//using tensorflow::tstring;
using tensorflow::Tensor;
using tensorflow::SavedModelBundle;
using tensorflow::SessionOptions;
using tensorflow::RunOptions;

class heaModel{
	private:
		SavedModelBundle bundle;
		SessionOptions sessOpt;
		RunOptions runOpt;
	public:
		void LoadModel(string, SessionOptions);
		void Predict(Tensor input, vector<Tensor>&);
};



void heaModel::LoadModel(string model_path, SessionOptions options){
	sessOpt = options;
	sessOpt.config.mutable_gpu_options()->set_allow_growth(true);
	auto status = tensorflow::LoadSavedModel(sessOpt, runOpt, model_path, {"serve"}, &bundle);
	if (!status.ok()){
		cout << "Error in LoadModel" << endl;
	}
}

void heaModel::Predict(Tensor input, vector<Tensor> &pred){
	const string input_node = "input_1:0";
	string output_node = "Identity:0";
//	const string input_node = "serving_default_input_1:0";
//	string output_node = "StatefulPartitionedCall:0";
	vector<std::pair<string, Tensor>> data = {{input_node, input}};
	auto status = this->bundle.GetSession()->Run(data, {output_node}, {}, &pred);
	if (!status.ok()){
		cout << "Error in Predict" << endl;
		cout << status.error_message();
	}

}
#endif
