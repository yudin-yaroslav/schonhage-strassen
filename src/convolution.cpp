#include "convolution.hpp"
#include "fft.hpp"

poly_ivector PWC_with_PROU(const poly_ivector &f, const poly_ivector &g, ivector w, uint64_t m) {
	// Data: Vectors f = [f_0, f_1, ..., f_{n-1}] and g = [g_0, g_1, ..., g_{n-1}]
	// with entries in R and an w in R that is a n-PROU.
	// Result: The positively wrapped convolution h+ = [h_0, ..., h_{n-1}] of f and g.

	if (f.size() != g.size()) {
		throw std::invalid_argument("PWC_with_PROU: Given `f` and `g` have different sizes");
	}
	const uint64_t n = f.size();

	poly_ivector a = FFT(f, w, m);
	poly_ivector b = FFT(g, w, m);

	poly_ivector c(n);
	for (size_t i = 0; i < n; i++) {
		c[i] = poly_mul_naive(a[i], b[i], m);
	}

	poly_ivector h = IFFT(c, w, m);

	return h;
}

poly_ivector NWC_with_PROU(const poly_ivector &f, const poly_ivector &g, ivector w, uint64_t m) {
	// Data: Vectors f = [f_0, f_1, ..., f_{n-1}] and g = [g_0, g_1, ..., g_{n-1}]
	// with entries in R and an w in R that is a 2n-PROU.
	// Result: The negatively wrapped convolution h- = [h_0, ..., h_{n-1}] of f and g.

	if (f.size() != g.size()) {
		throw std::invalid_argument("NMC_with_PROU: Given `f` and `g` have different sizes");
	}
	const uint64_t n = f.size();
	const ivector inv_w = poly_pow(w, n - 1, m);

	poly_ivector f_prime(n), g_prime(n);

	ivector pow_w = {1};
	for (size_t i = 0; i < n; ++i) {
		f_prime[i] = poly_mul_naive(pow_w, f[i], m);
		g_prime[i] = poly_mul_naive(pow_w, g[i], m);
		pow_w = poly_mul_naive(pow_w, w, m);
	}

	const ivector w_squared = poly_mul_naive(w, w, m);
	poly_ivector h_prime = PWC_with_PROU(f_prime, g_prime, w_squared, m);

	poly_ivector h(n);
	for (size_t i = 0; i < n; i++) {
		h[i] = poly_mul_naive(inv_w, h_prime[i], m);
	}

	return h;
}

ivector fast_NWC(const ivector &f, const ivector &g) {
	if (f.size() != g.size()) {
		throw std::invalid_argument("fast_NMC: Given `f` and `g` have different sizes");
	}
	const uint64_t n = f.size();
	uint64_t l = std::countr_zero(n);

	if ((n & (n - 1)) != 0) {
		throw std::invalid_argument("fast_NMC: Given `n` is not a power of 2");
	}

	const size_t k = (1ULL << ((l + 1) / 2));
	const size_t m = (1ULL << (l / 2));

	poly_ivector f_prime(k), g_prime(k);
	for (size_t i = 0; i < k; i++) {
		f_prime[i].resize(m);
		g_prime[i].resize(m);

		for (size_t j = 0; j < m; j++) {
			f_prime[i][j] = f[i * m + j];
			g_prime[i][j] = g[i * m + j];
		}
	}

	ivector w;
	if ((l & 1) == 0) {
		w = {0, 0, 1};
	} else {
		w = {0, 1};
	}

	poly_ivector h_prime = NWC_with_PROU(f_prime, g_prime, w, m);
	ivector h(n);

	// Evaluation of h(x, x^m)
	for (size_t i = 0; i < k; ++i) {
		const ivector &poly = h_prime[i];

		ivector tmp(n, 0);
		size_t dest = (i * m) % n;

		for (size_t j = 0; j < poly.size(); ++j) {
			tmp[dest + j] = poly[j];
		}

		if (i < k) {
			h = poly_add(h, tmp);
		} else {
			h = poly_sub(h, tmp);
		}
	}

	return h;
}
