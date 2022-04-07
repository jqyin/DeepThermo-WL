#ifndef SmartRedis_Client_H
#define SmartRedis_Client_H
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "client.h"


std::string encoder_name, decoder_name;
SmartRedis::Client* SRclient;

extern MPIState mpiState; 

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
		hea += std::string(alloyState.element[i+1]);
         
	encoder_name = "encoder_"+hea;
        std::string model = model_path + "/" + encoder_name + ".pb";
        (*SRclient).set_model_from_file(encoder_name+std::to_string(mpiState.myrank), model, "TF", "GPU", 
				1, 1, "vae", {"input_1"}, {"Identity"});
	std::cout << encoder_name << " loaded: " << (*SRclient).model_exists(encoder_name+std::to_string(mpiState.myrank)) << std::endl;

	decoder_name = "decoder_"+hea;
        model = model_path + "/" + decoder_name + ".pb";
        (*SRclient).set_model_from_file(decoder_name+std::to_string(mpiState.myrank), model, "TF", "GPU", 
				1, 1, "vae", {"input_2"}, {"Identity"});
	std::cout << decoder_name << " loaded: " << (*SRclient).model_exists(decoder_name+std::to_string(mpiState.myrank)) << std::endl;
}

#endif
