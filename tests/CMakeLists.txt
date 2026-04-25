enable_testing()

# Unit tests do not require MPI to be initialised; they exercise the pure
# state / parser / energy / RNG paths only. The deepthermo_core library is
# linked because the test TU references symbols defined there, but no MPI
# collectives are invoked.

add_executable(deepthermo_unit_tests
    unit/test_main.cc
    unit/test_parameter.cc
    unit/test_sim_context.cc
    unit/test_alloy_energy.cc
    unit/test_rand.cc
)
target_link_libraries(deepthermo_unit_tests PRIVATE deepthermo_core)
target_include_directories(deepthermo_unit_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/third_party
    ${CMAKE_SOURCE_DIR}/src
)
target_compile_definitions(deepthermo_unit_tests PRIVATE
    DEEPTHERMO_TESTS_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures"
)

add_test(NAME deepthermo_unit COMMAND deepthermo_unit_tests)

# Smoke test: shell-driven, runs hea-wl on a tiny config + freshly exported
# torch model and checks the expected output files appear. Skipped if
# python/torch aren't on PATH.
add_test(NAME deepthermo_smoke
    COMMAND ${CMAKE_COMMAND} -E env
        DEEPTHERMO_BUILD_DIR=${CMAKE_BINARY_DIR}
        DEEPTHERMO_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        bash ${CMAKE_CURRENT_SOURCE_DIR}/smoke/run.sh)
set_tests_properties(deepthermo_smoke PROPERTIES TIMEOUT 120)
