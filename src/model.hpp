#ifndef  MODEL_H
#define MODEL_H
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <tensorflow/core/public/session.h>
#include <tensorflow/core/protobuf/meta_graph.pb.h>
#include <tensorflow/core/public/session_options.h>
#include <tensorflow/cc/saved_model/loader.h>

using namespace std;
using tensorflow::Tensor;
using tensorflow::SavedModelBundle;
using tensorflow::SessionOptions;
using tensorflow::RunOptions;

enum PredictMode {encoder, decoder}; 

class heaModel{
	private:
		SavedModelBundle bundle;
		SessionOptions sessOpt;
		RunOptions runOpt;
	public:
		void LoadModel(string, SessionOptions);
		void Predict(Tensor input, vector<Tensor>&, PredictMode);
};



void heaModel::LoadModel(string model_path, SessionOptions options){
	sessOpt = options;
	sessOpt.config.mutable_gpu_options()->set_allow_growth(true);
	auto status = tensorflow::LoadSavedModel(sessOpt, runOpt, model_path, {"serve"}, &bundle);
	if (!status.ok()){
		cout << "Error in loading model" << endl;
	}
}

void heaModel::Predict(Tensor input, vector<Tensor> &pred, PredictMode mode){
        string input_node, output_node; 
        if (mode == encoder){ 	
		input_node = "input_1";
	}else if (mode == decoder){
		input_node = "input_2";
	}
	output_node = "Identity";
	vector<std::pair<string, Tensor>> data = {{input_node, input}};
	auto status = this->bundle.GetSession()->Run(data, {output_node}, {}, &pred);
	if (!status.ok()){
		cout << "Error in predict" << endl;
		cout << status.error_message();
	}

}
#endif
