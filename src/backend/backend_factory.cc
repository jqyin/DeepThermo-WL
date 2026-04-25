#include "backend/inference_backend.hpp"

#include <stdexcept>

#if defined(DEEPTHERMO_BACKEND_TORCH)
#include "backend/torch_backend.hpp"
#elif defined(DEEPTHERMO_BACKEND_TF)
#include "backend/tf_backend.hpp"
#elif defined(DEEPTHERMO_BACKEND_REDIS)
#include "backend/smartredis_backend.hpp"
#else
#error "No DEEPTHERMO_BACKEND_* macro defined; check CMake configure step."
#endif

std::unique_ptr<InferenceBackend> make_backend() {
#if defined(DEEPTHERMO_BACKEND_TORCH)
    return std::make_unique<TorchBackend>();
#elif defined(DEEPTHERMO_BACKEND_TF)
    return std::make_unique<TFBackend>();
#elif defined(DEEPTHERMO_BACKEND_REDIS)
    return std::make_unique<SmartRedisBackend>();
#endif
}
