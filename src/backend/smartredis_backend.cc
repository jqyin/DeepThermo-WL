#include "backend/smartredis_backend.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void SmartRedisBackend::load(const std::string& dir, int local_gpu, int rank,
                             const deepthermo::AlloyState& alloy) {
    // Co-located redis-server contract: each rank talks to a server on its
    // local node, port 6781 + (local_rank % gpus_per_node). The launcher
    // (e.g. examples/perlmutter/hea-wl.sb) starts these servers.
    std::string redis_host = "127.0.0.1:" + std::to_string(6781 + local_gpu);
    setenv("SSDB", redis_host.c_str(), 1);
    std::cout << "SmartRedisBackend: rank=" << rank << " SSDB=" << redis_host << std::endl;

    client_ = std::make_unique<SmartRedis::Client>(false);
    rank_str_ = std::to_string(rank);

    std::string hea;
    for (int i = 0; i < alloy.NE; ++i) hea += alloy.element[i + 1];

    encoder_name_ = "encoder_" + hea;
    {
        const std::string model = dir + "/" + encoder_name_ + ".pb";
        client_->set_model_from_file(encoder_name_ + rank_str_, model,
                                     "TF", "GPU", 1, 1, "vae",
                                     {"input_1"}, {"Identity"});
        if (!client_->model_exists(encoder_name_ + rank_str_)) {
            throw std::runtime_error("SmartRedisBackend: failed to load " + model);
        }
    }

    decoder_name_ = "decoder_" + hea;
    {
        const std::string model = dir + "/" + decoder_name_ + ".pb";
        client_->set_model_from_file(decoder_name_ + rank_str_, model,
                                     "TF", "GPU", 1, 1, "vae",
                                     {"input_2"}, {"Identity"});
        if (!client_->model_exists(decoder_name_ + rank_str_)) {
            throw std::runtime_error("SmartRedisBackend: failed to load " + model);
        }
    }
}

void SmartRedisBackend::encode(const float* input, int D, int NE, float* z) {
    const std::string z_key = "output_z_" + rank_str_;
    const std::string config_key = "input_config_" + rank_str_;
    client_->put_tensor(config_key, const_cast<float*>(input),
                        {1, static_cast<std::size_t>(D),
                         static_cast<std::size_t>(D),
                         static_cast<std::size_t>(D),
                         static_cast<std::size_t>(NE)},
                        SmartRedis::TensorType::flt,
                        SmartRedis::MemoryLayout::contiguous);
    client_->run_model(encoder_name_ + rank_str_, {config_key}, {z_key});
    client_->unpack_tensor(z_key, z, {3}, SmartRedis::TensorType::flt,
                           SmartRedis::MemoryLayout::contiguous);
}

void SmartRedisBackend::decode(const float* z, int D, int NE, float* output) {
    const std::string z_key = "input_z_" + rank_str_;
    const std::string config_key = "output_config_" + rank_str_;
    client_->put_tensor(z_key, const_cast<float*>(z), {1, 3},
                        SmartRedis::TensorType::flt,
                        SmartRedis::MemoryLayout::contiguous);
    client_->run_model(decoder_name_ + rank_str_, {z_key}, {config_key});
    client_->unpack_tensor(config_key, output,
                           {static_cast<std::size_t>(D) * D * D * NE},
                           SmartRedis::TensorType::flt,
                           SmartRedis::MemoryLayout::contiguous);
}
