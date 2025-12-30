#pragma once

#include "utils.hpp"

// Given two vectors f = [f_0, f_1, ..., f_{n-1}] and g = [g_0, g_1, ..., g_{n-1}]
// with each f_i, g_i in R, the positively wrapped convolution(PWC) of f and g
// is given by the vector h+ = [h+_0, ..., h+_{n-1}] defined as:
// h+_l = sum_{i=0}^{n-1} (f_i g_{l-i} + f_i g_{n+l-i}).
// Similarly, the negatively wrapped convolution(NWC) of f and g is given
// by the vector h- = [h-_0, ..., h-_{n-1}] defined as
// h-_l = sum_{i=0}^{n-1} (f_i g_{l-i} - f_i g_{n+l - i}).
//
// Note that
// h+(x) = f(x) g(x) mod (x^n - 1)
// h-(x) = f(x) g(x) mod (x^n + 1)
//
// WARNING: the functions PWC_with_PROU and NWC_with_PROU
// are written to work in the ring of polynomials
// Z[x]/(x^{2m} + 1), not in the field of residues modulo some prime.

poly_ivector PWC_with_PROU(const poly_ivector &f, const poly_ivector &g, ivector w, uint64_t m);
poly_ivector NWC_with_PROU(const poly_ivector &f, const poly_ivector &g, ivector w, uint64_t m);
ivector naive_NWC(const ivector &f, const ivector &g);
ivector fast_NWC(const ivector &f, const ivector &g);
ivector full_product(const ivector &f, const ivector &g);
