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

        # Some ROCm LibTorch wheels bake the *unversioned* /opt/rocm/lib path
        # into c10_hip's/torch_hip's INTERFACE_LINK_LIBRARIES (from whatever
        # build host produced the wheel), instead of the ROCM_PATH that was
        # actually resolved here. HPC systems that keep only versioned
        # /opt/rocm-<ver> trees (no /opt/rocm symlink) then fail at build.make
        # generation time with "No rule to make target /opt/rocm/lib/...".
        # Rewrite those entries to the ROCm install CMake actually found.
        if((TARGET c10_hip OR TARGET torch_hip) AND NOT EXISTS "/opt/rocm/lib")
            if(DEFINED ENV{ROCM_PATH})
                set(_deepthermo_rocm_lib_dir "$ENV{ROCM_PATH}/lib")
            elseif(HIP_ROOT_DIR)
                set(_deepthermo_rocm_lib_dir "${HIP_ROOT_DIR}/lib")
            endif()
            if(_deepthermo_rocm_lib_dir)
                foreach(_deepthermo_hip_tgt c10_hip torch_hip)
                    if(TARGET ${_deepthermo_hip_tgt})
                        get_target_property(_deepthermo_hip_libs ${_deepthermo_hip_tgt} INTERFACE_LINK_LIBRARIES)
                        if(_deepthermo_hip_libs)
                            string(REPLACE "/opt/rocm/lib/" "${_deepthermo_rocm_lib_dir}/"
                                   _deepthermo_hip_libs "${_deepthermo_hip_libs}")
                            set_target_properties(${_deepthermo_hip_tgt} PROPERTIES
                                INTERFACE_LINK_LIBRARIES "${_deepthermo_hip_libs}")
                        endif()
                    endif()
                endforeach()
            endif()
            unset(_deepthermo_rocm_lib_dir)
        endif()

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
