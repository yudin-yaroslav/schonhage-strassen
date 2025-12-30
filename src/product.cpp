#include "product.hpp"
#include "fft.hpp"

poly_ivector PWC_with_PROU(const poly_ivector &f, const poly_ivector &g, ivector w, uint64_t m) {
	// Data: Vectors f = [f_0, f_1, ..., f_{n-1}] and g = [g_0, g_1, ..., g_{n-1}]
	// with entries in R and an w in R that is a n-PROU.
	// Result: The positively wrapped convolution h+ = [h_0, ..., h_{n-1}] of f and g.

	if (f.size() != g.size()) {
		throw std::invalid_argument("PWC_with_PROU: Given `f` and `g` have different sizes");
	}
	const size_t n = f.size();

	poly_ivector a = FFT(f, w, m);
	poly_ivector b = FFT(g, w, m);

	poly_ivector c(n);
	for (size_t i = 0; i < n; i++) {
		ivector a_padded = a[i];
		ivector b_padded = b[i];
		a_padded.resize(2 * m, 0);
		b_padded.resize(2 * m, 0);

		c[i] = fast_NWC(a_padded, b_padded);
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
	const size_t n = f.size();

	uint64_t e;
	if (n == m) {
		e = 2;
	} else if (n == 2 * m) {
		e = 1;
	} else {
		throw std::invalid_argument("NWC_with_PROU: Invalid size relationship between n and m");
	}

	poly_ivector f_prime(n), g_prime(n);
	for (size_t i = 0; i < n; ++i) {
		f_prime[i] = poly_mul_by_x_power(f[i], i * e, m);
		g_prime[i] = poly_mul_by_x_power(g[i], i * e, m);
	}

	ivector w_squared = poly_mul_by_x_power(w, e, m);
	poly_ivector h_prime = PWC_with_PROU(f_prime, g_prime, w_squared, m);
	poly_ivector h(n);
	for (size_t i = 0; i < n; ++i) {
		h[i] = poly_mul_by_x_power(h_prime[i], (2 * n - 1) * i * e, m);
	}
	return h;
}

ivector naive_NWC(const ivector &f, const ivector &g) {
	const size_t n = f.size();
	ivector full_conv(2 * n - 1, 0);

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			full_conv[i + j] += f[i] * g[j];
		}
	}

	ivector result(n, 0);
	for (size_t i = 0; i < n; i++) {
		result[i] = full_conv[i];
	}
	for (size_t i = n; i < 2 * n - 1; i++) {
		result[i - n] -= full_conv[i];
	}

	return result;
}

ivector fast_NWC(const ivector &f, const ivector &g) {
	// Data: Vectors f = [f_0, f_1, ..., f_{n-1}] and g = [g_0, g_1, ..., g_{n-1}]
	// with entries in R and n = 2^l.
	// Result: The negatively wrapped convolution h- = [h_0, ..., h_{n-1}] of f and g.

	if (f.size() != g.size()) {
		throw std::invalid_argument("fast_NWC: Given `f` and `g` have different sizes");
	}

	const size_t n = f.size();

	if ((n & (n - 1)) != 0) {
		throw std::invalid_argument("fast_NWC: Given `n` is not a power of 2");
	}

	// Base case
	if (n <= 64) {
		return naive_NWC(f, g);
	}

	const size_t l = std::countr_zero(n);
	const size_t m = (1ULL << (l / 2));
	const size_t k = (1ULL << ((l + 1) / 2));

	poly_ivector f_tilde(k), g_tilde(k);
	for (size_t i = 0; i < k; i++) {
		f_tilde[i].resize(m, 0);
		g_tilde[i].resize(m, 0);

		for (size_t j = 0; j < m; j++) {
			f_tilde[i][j] = f[i * m + j];
			g_tilde[i][j] = g[i * m + j];
		}
	}

	for (size_t i = 0; i < k; i++) {
		f_tilde[i].resize(2 * m, 0);
		g_tilde[i].resize(2 * m, 0);
	}

	ivector w;
	uint64_t phi_deg = 2 * m;

	if (l % 2 == 0) {
		// l even: k = m, w = x^2 (2k-PROU = 2m-PROU)
		w = monomial_from_exponent(2, phi_deg);
	} else {
		// l odd: k = 2m, w = x (2k-PROU = 4m-PROU)
		w = monomial_from_exponent(1, phi_deg);
	}

	poly_ivector h_tilde = NWC_with_PROU(f_tilde, g_tilde, w, m);

	ivector h(n, 0);

	for (size_t i = 0; i < k; i++) {
		const ivector &poly = h_tilde[i];

		for (size_t j = 0; j < 2 * m; j++) {
			int64_t coeff = poly[j];
			if (coeff == 0)
				continue;

			int64_t exp = static_cast<int64_t>(i) * m + j;

			while (exp >= static_cast<int64_t>(n)) {
				exp -= n;
				coeff = -coeff;
			}

			h[exp] += coeff;
		}
	}

	return h;
}

ivector full_product(const ivector &f, const ivector &g) {
	const size_t n = f.size();
	const size_t padded_len = 2 * n;

	ivector f_padded(padded_len, 0);
	ivector g_padded(padded_len, 0);

	for (size_t i = 0; i < n; i++) {
		f_padded[i] = f[i];
		g_padded[i] = g[i];
	}

	ivector result = fast_NWC(f_padded, g_padded);

	result.resize(2 * n - 1, 0);
	return result;
}
