cmake_minimum_required(VERSION 3.18)
project(DeepThermo
    VERSION 2.0.0
    DESCRIPTION "Deep learning-accelerated parallel Monte Carlo sampling for high-entropy alloys"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---- options --------------------------------------------------------------

set(DEEPTHERMO_BACKEND "torch" CACHE STRING
    "Inference backend: torch (default, LibTorch), tf (TensorFlow C++), redis (SmartRedis).")
set_property(CACHE DEEPTHERMO_BACKEND PROPERTY STRINGS torch tf redis)

set(DEEPTHERMO_PLATFORM "generic" CACHE STRING
    "Target platform: generic (default), frontier (OLCF, MI250X), perlmutter (NERSC, A100).")
set_property(CACHE DEEPTHERMO_PLATFORM PROPERTY STRINGS generic frontier perlmutter)

option(DEEPTHERMO_BUILD_TESTS "Build unit + smoke tests" ON)
option(DEEPTHERMO_ENABLE_GLOBAL_UPDATE "Enable VAE-driven global update (-DGLOBAL_UPDATE)" ON)
option(DEEPTHERMO_ENABLE_RING_ALLREDUCE "Use ring allreduce instead of MPI_Allreduce" OFF)

message(STATUS "DeepThermo backend  : ${DEEPTHERMO_BACKEND}")
message(STATUS "DeepThermo platform : ${DEEPTHERMO_PLATFORM}")

# ---- dependencies ---------------------------------------------------------

find_package(MPI REQUIRED COMPONENTS CXX)

add_library(deepthermo_third_party INTERFACE)
target_include_directories(deepthermo_third_party
    INTERFACE ${CMAKE_SOURCE_DIR}/third_party)

include(${CMAKE_SOURCE_DIR}/cmake/backends.cmake)
deepthermo_configure_backend(${DEEPTHERMO_BACKEND})

# ---- platform-specific tuning --------------------------------------------

add_library(deepthermo_platform INTERFACE)
if(DEEPTHERMO_PLATFORM STREQUAL "frontier")
    target_compile_definitions(deepthermo_platform INTERFACE DEEPTHERMO_PLATFORM_FRONTIER)
elseif(DEEPTHERMO_PLATFORM STREQUAL "perlmutter")
    target_compile_definitions(deepthermo_platform INTERFACE DEEPTHERMO_PLATFORM_PERLMUTTER)
endif()

# ---- targets --------------------------------------------------------------

add_subdirectory(src)

if(DEEPTHERMO_BUILD_TESTS)
    enable_testing()
    if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
        add_subdirectory(tests)
    endif()
endif()
