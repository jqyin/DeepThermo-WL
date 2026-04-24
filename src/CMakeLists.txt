set(DEEPTHERMO_SOURCES
    alloy.cc
    allreduce.cc
    main.cc
    parameter.cc
    pt.cc
    rand.cc
    wanglandau.cc
)

add_executable(hea-wl ${DEEPTHERMO_SOURCES})

target_link_libraries(hea-wl
    PRIVATE
        deepthermo_backend
        deepthermo_platform
        deepthermo_third_party
        MPI::MPI_CXX
)

if(DEEPTHERMO_ENABLE_GLOBAL_UPDATE)
    target_compile_definitions(hea-wl PRIVATE GLOBAL_UPDATE)
endif()

if(DEEPTHERMO_ENABLE_RING_ALLREDUCE)
    target_compile_definitions(hea-wl PRIVATE COMM_RING)
endif()

# Legacy synonym while the refactor is in flight; src/main.hpp still checks
# these names. They go away in the cleanup pass (Task 7 / Task 9).
if(DEEPTHERMO_BACKEND STREQUAL "tf")
    target_compile_definitions(hea-wl PRIVATE TF_BACKEND)
endif()
if(DEEPTHERMO_PLATFORM STREQUAL "frontier")
    target_compile_definitions(hea-wl PRIVATE CRUSHER)
endif()

target_compile_options(hea-wl PRIVATE
    -Wno-narrowing
    -rdynamic
)

set_target_properties(hea-wl PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
