#pragma once

#include "utils.hpp"

ivector NTT(const ivector &coeffs, uint8_t k, uint64_t mod, uint64_t w_n);
ivector INTT(const ivector &values, uint8_t k, uint64_t mod, uint64_t w_n);
