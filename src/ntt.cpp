#include "ntt.hpp"
#include <cmath>

ivector NTT(const ivector &coeffs, uint8_t k, uint64_t mod, uint64_t w_n) {
	// ord_mod(w_n) = n = 2^k
	// mod is a prime
	// mod - 1 is divisible by n

	uint64_t n = (1ULL << k);

	if (k == 0) {
		return {coeffs[0] % mod};
	}

	if (coeffs.size() != n) {
		throw std::invalid_argument("Given polynomial coeffs' number is not equal to n");
	}
	if (pow_mod(w_n, n, mod) != 1) {
		throw std::invalid_argument("Given w_n is not an n-th root of unity modulo mod");
	}

	ivector coeffs_even(n / 2), coeffs_odd(n / 2);
	for (size_t i = 0; i < n / 2; ++i) {
		coeffs_even[i] = coeffs[2 * i];
		coeffs_odd[i] = coeffs[2 * i + 1];
	}

	uint64_t w_squared = mul_mod(w_n, w_n, mod);
	ivector values_even = NTT(coeffs_even, k - 1, mod, w_squared);
	ivector values_odd = NTT(coeffs_odd, k - 1, mod, w_squared);

	ivector values(n);

	uint64_t w = 1;
	for (size_t i = 0; i < n / 2; ++i) {
		uint64_t temp = mul_mod(w, values_odd[i], mod);

		values[i] = add_mod(values_even[i], temp, mod);
		values[i + n / 2] = sub_mod(values_even[i], temp, mod);

		w = mul_mod(w, w_n, mod);
	}

	return values;
}

ivector INTT(const ivector &values, uint8_t k, uint64_t mod, uint64_t w_n) {
	uint64_t n = (1ULL << k);
	uint64_t inv_n = pow_mod(n, mod - 2, mod);
	uint64_t inv_w = pow_mod(w_n, mod - 2, mod);

	ivector coeffs_temp = NTT(values, k, mod, inv_w);

	ivector coeffs(n);
	for (size_t i = 0; i < n; i++) {
		coeffs[i] = mul_mod(inv_n, coeffs_temp[i], mod);
	}
	return coeffs;
}
