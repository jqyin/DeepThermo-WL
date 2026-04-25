/**
 * @file main.hpp
 * @brief Legacy shim — pulls in sim_context.hpp.
 *
 * Historically this header declared the four extern state singletons
 * plus a pile of `#define`'d alloy constants. All of that now lives in
 * sim_context.hpp as proper classes with runtime fields. New code should
 * include sim_context.hpp directly; this header is retained only because
 * existing translation units reference `#include "main.hpp"`.
 */

#pragma once

#include "sim_context.hpp"
