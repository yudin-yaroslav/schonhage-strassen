#include "fft.hpp"

// WARNING: R = Z[x]/(x^{2m} + 1)

poly_ivector FFT(const poly_ivector &f, ivector w, uint64_t m) {
	// Data: f_0, f_1, ..., f_{n-1} in R, and w in R that is an n-PROU
	// Result: (f(1), f(w), . . . , f(w^{n−1}))

	const uint64_t n = f.size();

	if (n == 1) {
		return {f[0]};
	}

	if ((n & (n - 1)) != 0) {
		throw std::invalid_argument("FFT: Given `n` is not a power of 2");
	}

	poly_ivector f_even(n / 2), f_odd(n / 2);
	for (size_t i = 0; i < n / 2; ++i) {
		f_even[i] = f[2 * i];
		f_odd[i] = f[2 * i + 1];
	}

	const ivector w_squared = poly_mul_naive(w, w, m);
	poly_ivector values_even = FFT(f_even, w_squared, m);
	poly_ivector values_odd = FFT(f_odd, w_squared, m);

	poly_ivector values(n);

	ivector pow_w(2 * m, 0);
	pow_w[0] = 1;

	for (size_t i = 0; i < n / 2; ++i) {
		ivector temp = poly_mul_naive(pow_w, values_odd[i], m);

		values[i] = poly_add(values_even[i], temp);
		values[i + n / 2] = poly_sub(values_even[i], temp);

		pow_w = poly_mul_naive(pow_w, w, m);
	}

	return values;
}

poly_ivector IFFT(const poly_ivector &f, ivector w, uint64_t m) {
	const uint64_t n = f.size();

	if ((n & (n - 1)) != 0) {
		throw std::invalid_argument("IFFT: Given `n` is not a power of 2");
	}

	const ivector inv_w = poly_pow(w, n - 1, m);

	poly_ivector h = FFT(f, inv_w, m);

	for (size_t i = 0; i < h.size(); i++) {
		h[i] = poly_div_scalar(h[i], n);
	}
	return h;
}
