#include "convolution.hpp"
#include "ntt.hpp"

ivector PWC_with_PROU(const ivector &f, const ivector &g, uint64_t w, uint64_t mod) {
	// Data: Vectors f = [f_0, f_1, ..., f_{n-1}] and g = [g_0, g_1, ..., g_{n-1}]
	// with entries in R and an w in R that is a n-PROU.
	// Result: The positively wrapped convolution h+ = [h_0, ..., h_{n-1}] of f and g.

	if (f.size() != g.size()) {
		throw std::invalid_argument("Given `f` and `g` have different sizes");
	}
	const uint64_t n = f.size();

	ivector a = NTT(f, w, mod);
	ivector b = NTT(g, w, mod);

	ivector c(n);
	for (size_t i = 0; i < n; i++) {
		c[i] = mul_mod(a[i], b[i], mod);
	}

	ivector h = INTT(c, w, mod);

	return h;
}

ivector NWC_with_PROU(const ivector &f, const ivector &g, uint64_t w, uint64_t mod) {
	// Data: Vectors f = [f_0, f_1, ..., f_{n-1}] and g = [g_0, g_1, ..., g_{n-1}]
	// with entries in R and an w in R that is a 2n-PROU.
	// Result: The negatively wrapped convolution h- = [h_0, ..., h_{n-1}] of f and g.

	if (f.size() != g.size()) {
		throw std::invalid_argument("Given `f` and `g` have different sizes");
	}
	const uint64_t n = f.size();
	const uint64_t inv_w = pow_mod(w, mod - 2, mod);

	ivector f_prime(n), g_prime(n);
	uint64_t pow_w = 1;
	for (size_t i = 0; i < n; ++i) {
		f_prime[i] = mul_mod(pow_w, f[i], mod);
		g_prime[i] = mul_mod(pow_w, g[i], mod);
		pow_w = mul_mod(pow_w, w, mod);
	}

	const uint64_t w_squared = mul_mod(w, w, mod);
	ivector h_prime = PWC_with_PROU(f_prime, g_prime, w_squared, mod);

	ivector h(n);
	for (size_t i = 0; i < n; i++) {
		h[i] = mul_mod(inv_w, h_prime[i], mod);
	}

	return h;
}

ivector fast_NWC(const ivector &f, const ivector &g, uint64_t mod) {
	if (f.size() != g.size()) {
		throw std::invalid_argument("Given `f` and `g` have different sizes");
	}
	const uint64_t n = f.size();
	uint64_t l = std::countr_zero(n);

	if ((n & (n - 1)) != 0) {
		throw std::invalid_argument("Given `n` is not a power of 2");
	}

	const uint64_t k = (1ULL << ((l + 1) / 2));
	const uint64_t m = (1ULL << (l / 2));

	uint64_t w;
	if ((l & 1) == 0) {
		w = x
	}
}
