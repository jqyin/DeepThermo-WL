// main.hpp — legacy shim header.
//
// Historically this header declared the four extern state singletons plus a
// pile of #define'd alloy constants. All of that now lives in sim_context.hpp
// as proper classes with runtime fields. New code should include
// sim_context.hpp directly; this header is kept because existing translation
// units all #include "main.hpp" and we change one concern at a time.

#pragma once

#include "sim_context.hpp"
