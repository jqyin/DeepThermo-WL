#ifndef SmartRedis_Client_H
#define SmartRedis_Client_H
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "client.h"

#define GPUperNode 6

std::string model_name;
SmartRedis::Client* SRclient;

extern int myrank; 
extern int nprocs; 

void LoadClientModel(std::string model_path){
        
    	const char* rank = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    	int local_rank = atoi(rank);

	std::string redis_host = "127.0.0.1:" + std::to_string(local_rank%GPUperNode+6781);
	setenv("SSDB", redis_host.c_str(), 1);

        const char* SSDB = getenv("SSDB");
        std::cout << "local_rank " << local_rank << " SSDB: " << SSDB << std::endl;

	SRclient = new SmartRedis::Client(false);

        std::string hea=""; 
 	for(int i=0; i< NE; i++)
		hea += std::string(element[i+1]);
         
	model_name = "vae_"+hea;
        std::string model = model_path + "/" + model_name + ".pb";
        (*SRclient).set_model_from_file(model_name+std::to_string(myrank), model, "TF", "GPU", 
				1, SH*NE*(NE-1)/2, "vae", {"input_2"}, {"conv3d_transpose_4_1/Sigmoid"});
	std::cout << model_name << " loaded: " << (*SRclient).model_exists(model_name+std::to_string(myrank)) << std::endl;
}

#endif
