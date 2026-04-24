// client.hpp — SmartRedis inference-server client (legacy backend).
//
// Only included when DEEPTHERMO_BACKEND=redis. Inline globals (C++17)
// eliminate the ODR issue the pre-refactor code papered over with
// --allow-multiple-definition. Task 4 moves this into a proper
// SmartRedisBackend implementing the InferenceBackend interface.

#pragma once

#include <cstdlib>
#include <iostream>
#include <string>

#include "client.h"

#include "sim_context.hpp"

inline std::string encoder_name;
inline std::string decoder_name;
inline SmartRedis::Client* SRclient = nullptr;

inline void LoadClientModel(std::string model_path) {
    const char* rank = std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    const int local_rank = rank ? std::atoi(rank) : 0;

    std::string redis_host = "127.0.0.1:" +
        std::to_string(local_rank % mpiState.gpus_per_node + 6781);
    setenv("SSDB", redis_host.c_str(), 1);

    const char* SSDB = std::getenv("SSDB");
    std::cout << "local_rank " << local_rank << " SSDB: " << SSDB << std::endl;

    SRclient = new SmartRedis::Client(false);

    std::string hea;
    for (int i = 0; i < alloyState.NE; ++i) {
        hea += alloyState.element[i + 1];
    }

    encoder_name = "encoder_" + hea;
    std::string model = model_path + "/" + encoder_name + ".pb";
    (*SRclient).set_model_from_file(encoder_name + std::to_string(mpiState.myrank), model,
                                    "TF", "GPU", 1, 1, "vae", {"input_1"}, {"Identity"});
    std::cout << encoder_name << " loaded: "
              << (*SRclient).model_exists(encoder_name + std::to_string(mpiState.myrank))
              << std::endl;

    decoder_name = "decoder_" + hea;
    model = model_path + "/" + decoder_name + ".pb";
    (*SRclient).set_model_from_file(decoder_name + std::to_string(mpiState.myrank), model,
                                    "TF", "GPU", 1, 1, "vae", {"input_2"}, {"Identity"});
    std::cout << decoder_name << " loaded: "
              << (*SRclient).model_exists(decoder_name + std::to_string(mpiState.myrank))
              << std::endl;
}
