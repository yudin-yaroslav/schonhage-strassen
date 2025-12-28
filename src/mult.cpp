#include "mult.hpp"
#include "ntt.hpp"
#include <bit>

ivector multiply(ivector a, ivector b) {
	// n = 2^k
	// ord_p(g) = 2n = 2^{k+1}

	if ((a.size() & (a.size() - 1)) != 0 || (b.size() & (b.size() - 1)) != 0) {
		throw "The size of input number is not a power of 2";
	}
	uint64_t max_degree = a.size() + b.size() + 1;

	uint8_t k = std::bit_width(max_degree - 1);
	uint64_t n = (1ULL << k);

	a.resize(2 * n);
	b.resize(2 * n);

	uint64_t prime = ntt_primes[k + 1];
	uint64_t order_2n = ntt_order_k[k + 1];

	ivector a_values = NTT(a, k + 1, prime, order_2n);
	ivector b_values = NTT(b, k + 1, prime, order_2n);

	ivector a_coeffs_cycle = INTT(a_values, k + 1, prime, order_2n);

	for (int i = 0; i < 10; i++) {
		std::cout << a[i] % prime << " " << a_coeffs_cycle[i] << std::endl;
	}
	// for (int i = 0; i < max_degree; i++) {
	// 	ab_values[i] = a_values[i] * b_values[i];
	// }

	// ivector ab = IFFT(&ab_values, max_degree);

	ivector ab_values(2 * n);
	return ab_values;
}
