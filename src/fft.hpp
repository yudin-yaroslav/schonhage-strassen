#pragma once

#include "utils.hpp"

// WARNING: R = Z[x]/(x^{2m} + 1)

// Number Theoretic Transform
// Suppose w in R is an n-PROU in Z_m, where n is a power of 2.
// The Number Theoretic Transform (NTT) with respect to w
// of a polynomial f(x) in R[x] with deg(f) < n is the tuple
// of evaluations:
// NTT_w(f) = (f(1), f(w), ..., f(w^{n-1}))

poly_ivector FFT(const poly_ivector &f, ivector w, uint64_t m);
poly_ivector IFFT(const poly_ivector &f, ivector w, uint64_t m);
