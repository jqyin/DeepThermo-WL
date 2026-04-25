set(DEEPTHERMO_CORE_SOURCES
    alloy.cc
    allreduce.cc
    parameter.cc
    pt.cc
    rand.cc
    sim_context.cc
    wanglandau.cc
    backend/backend_factory.cc
)

if(DEEPTHERMO_BACKEND STREQUAL "torch")
    list(APPEND DEEPTHERMO_CORE_SOURCES backend/torch_backend.cc)
elseif(DEEPTHERMO_BACKEND STREQUAL "tf")
    list(APPEND DEEPTHERMO_CORE_SOURCES backend/tf_backend.cc)
elseif(DEEPTHERMO_BACKEND STREQUAL "redis")
    list(APPEND DEEPTHERMO_CORE_SOURCES backend/smartredis_backend.cc)
endif()

add_library(deepthermo_core STATIC ${DEEPTHERMO_CORE_SOURCES})

target_include_directories(deepthermo_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(deepthermo_core
    PUBLIC
        deepthermo_backend
        deepthermo_platform
        deepthermo_third_party
        MPI::MPI_CXX
)
target_compile_options(deepthermo_core PRIVATE -Wno-narrowing)
target_link_options(deepthermo_core PUBLIC -rdynamic)

if(DEEPTHERMO_ENABLE_GLOBAL_UPDATE)
    target_compile_definitions(deepthermo_core PUBLIC GLOBAL_UPDATE)
endif()
if(DEEPTHERMO_ENABLE_RING_ALLREDUCE)
    target_compile_definitions(deepthermo_core PUBLIC COMM_RING)
endif()

add_executable(hea-wl main.cc)
target_link_libraries(hea-wl PRIVATE deepthermo_core)

set_target_properties(hea-wl PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
