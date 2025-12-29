#include "ntt.hpp"

ivector NTT(const ivector &f, uint64_t w, uint64_t mod) {
	// Data: f_0, f_1, ..., f_{n-1} in R, and w in R that is an n-PROU
	// Result: (f(1), f(w), . . . , f(w^{n−1}))

	const uint64_t n = f.size();

	if (n == 1) {
		return {f[0]};
	}

	if ((n & (n - 1)) != 0) {
		throw std::invalid_argument("Given `n` is not a power of 2");
	}
	if (!is_prou(w, n, mod)) {
		throw std::invalid_argument("Given `w` is not an n-PROU modulo `mod`");
	}

	ivector f_even(n / 2), f_odd(n / 2);
	for (size_t i = 0; i < n / 2; ++i) {
		f_even[i] = f[2 * i];
		f_odd[i] = f[2 * i + 1];
	}

	const uint64_t w_squared = mul_mod(w, w, mod);
	ivector values_even = NTT(f_even, w_squared, mod);
	ivector values_odd = NTT(f_odd, w_squared, mod);

	ivector values(n);

	uint64_t pow_w = 1;
	for (size_t i = 0; i < n / 2; ++i) {
		uint64_t temp = mul_mod(pow_w, values_odd[i], mod);

		values[i] = add_mod(values_even[i], temp, mod);
		values[i + n / 2] = sub_mod(values_even[i], temp, mod);

		w = mul_mod(pow_w, w, mod);
	}

	return values;
}

ivector INTT(const ivector &f, uint64_t w, uint64_t mod) {
	const uint64_t n = f.size();

	const uint64_t inv_n = pow_mod(n, mod - 2, mod);
	const uint64_t inv_w = pow_mod(w, mod - 2, mod);

	ivector coeffs_temp = NTT(f, inv_w, mod);

	ivector coeffs(n);
	for (size_t i = 0; i < n; i++) {
		coeffs[i] = mul_mod(inv_n, coeffs_temp[i], mod);
	}
	return coeffs;
}
