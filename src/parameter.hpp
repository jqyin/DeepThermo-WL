#pragma once

// Loads <filename> (a TOML configuration) into the global SimContext fields.
// Throws std::runtime_error on parse / validation failure.
void ReadInput(const char* filename);
