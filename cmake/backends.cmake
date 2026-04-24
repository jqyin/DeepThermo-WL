# deepthermo_configure_backend(<backend>)
#
# Defines an INTERFACE target `deepthermo_backend` that, when linked, pulls in
# the headers + libraries needed to build against the selected inference
# backend. Callers set definitions via target_compile_definitions on this
# target so downstream code can branch at compile time.
#
# Exposes:
#   DEEPTHERMO_BACKEND_DEFINE   - the DEEPTHERMO_BACKEND_<X> preprocessor define
#   deepthermo_backend          - INTERFACE target to link

function(deepthermo_configure_backend backend)
    add_library(deepthermo_backend INTERFACE)

    if(backend STREQUAL "torch")
        find_package(Torch REQUIRED)
        # Torch's CMake config sets TORCH_CXX_FLAGS which may include
        # -D_GLIBCXX_USE_CXX11_ABI=<0|1>; propagate as INTERFACE.
        if(TORCH_CXX_FLAGS)
            separate_arguments(_torch_flags UNIX_COMMAND "${TORCH_CXX_FLAGS}")
            target_compile_options(deepthermo_backend INTERFACE ${_torch_flags})
        endif()
        target_link_libraries(deepthermo_backend INTERFACE ${TORCH_LIBRARIES})
        target_compile_definitions(deepthermo_backend INTERFACE DEEPTHERMO_BACKEND_TORCH)
        set(DEEPTHERMO_BACKEND_DEFINE "DEEPTHERMO_BACKEND_TORCH" PARENT_SCOPE)

    elseif(backend STREQUAL "tf")
        # TensorFlow C++ has no official CMake config; TF_DIR must point to a
        # TensorFlow install dir containing include/ and libtensorflow_cc.so.
        if(NOT DEFINED TF_DIR AND DEFINED ENV{TF_DIR})
            set(TF_DIR $ENV{TF_DIR})
        endif()
        if(NOT TF_DIR)
            message(FATAL_ERROR "BACKEND=tf requires TF_DIR to be set (via -DTF_DIR=... or env).")
        endif()
        target_include_directories(deepthermo_backend SYSTEM INTERFACE ${TF_DIR}/include)
        target_link_libraries(deepthermo_backend INTERFACE
            ${TF_DIR}/libtensorflow_cc.so
            ${TF_DIR}/libtensorflow_framework.so)
        target_link_options(deepthermo_backend INTERFACE
            "LINKER:-rpath,${TF_DIR}"
            "LINKER:--allow-multiple-definition")
        target_compile_definitions(deepthermo_backend INTERFACE
            DEEPTHERMO_BACKEND_TF
            _GLIBCXX_USE_CXX11_ABI=1)
        set(DEEPTHERMO_BACKEND_DEFINE "DEEPTHERMO_BACKEND_TF" PARENT_SCOPE)

    elseif(backend STREQUAL "redis")
        if(NOT DEFINED REDIS_DIR AND DEFINED ENV{REDIS_DIR})
            set(REDIS_DIR $ENV{REDIS_DIR})
        endif()
        if(NOT REDIS_DIR)
            message(FATAL_ERROR "BACKEND=redis requires REDIS_DIR to be set (via -DREDIS_DIR=... or env).")
        endif()
        target_include_directories(deepthermo_backend SYSTEM INTERFACE ${REDIS_DIR}/include)
        target_link_libraries(deepthermo_backend INTERFACE ${REDIS_DIR}/lib/libsmartredis.so)
        target_link_options(deepthermo_backend INTERFACE "LINKER:-rpath,${REDIS_DIR}/lib")
        target_compile_definitions(deepthermo_backend INTERFACE DEEPTHERMO_BACKEND_REDIS)
        set(DEEPTHERMO_BACKEND_DEFINE "DEEPTHERMO_BACKEND_REDIS" PARENT_SCOPE)

    else()
        message(FATAL_ERROR "Unknown backend '${backend}'. Choose torch|tf|redis.")
    endif()
endfunction()
