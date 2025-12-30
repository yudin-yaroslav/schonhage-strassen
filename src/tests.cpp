#include "fft.hpp"
#include "product.hpp"
#include "utils.hpp"
#include <format>
#include <random>

void test_poly_add() {
	const std::string name = "test_poly_add";

	{
		ivector f1 = {1, 2, 3};
		ivector g1 = {4, 5, 6};
		ivector expect1 = {5, 7, 9};
		ivector got1 = poly_add(f1, g1);
		if (!eq_ivector(got1, expect1)) {
			fail_expectation(name + " correctness", &got1, &expect1);
			return;
		}
	}

	{
		ivector f2 = {-5, 0, 10};
		ivector g2 = {5, 0, -10};
		ivector expect2 = {0, 0, 0};
		ivector got2 = poly_add(f2, g2);
		if (!eq_ivector(got2, expect2)) {
			fail_expectation(name + " neg-zero", nullptr, nullptr);
			return;
		}
	}

	try {
		ivector f3 = {1, 2};
		ivector g3 = {1, 2, 3};
		(void)poly_add(f3, g3);
		fail_expectation(name + " size-mismatch (no exception thrown)", nullptr, nullptr);
		return;
	} catch (const std::invalid_argument &) {
		// expected
	}

	pass(name);
}

void test_poly_sub() {
	const std::string name = "test_poly_sub";

	{
		ivector f1 = {5, 7, 9};
		ivector g1 = {1, 2, 3};
		ivector expect1 = {4, 5, 6};
		ivector got1 = poly_sub(f1, g1);
		if (!eq_ivector(got1, expect1)) {
			fail_expectation(name + " basic", &got1, &expect1);
			return;
		}
	}

	{
		ivector f2 = {0, 0, 0};
		ivector g2 = {1, -2, 3};
		ivector expect2 = {-1, 2, -3};
		ivector got2 = poly_sub(f2, g2);
		if (!eq_ivector(got2, expect2)) {
			fail_expectation(name + " negative", &got2, &expect2);
			return;
		}
	}

	try {
		ivector f3 = {1};
		ivector g3 = {1, 2};
		(void)poly_sub(f3, g3);
		fail_expectation(name + " size-mismatch (no exception thrown)", nullptr, nullptr);
		return;
	} catch (const std::invalid_argument &) {
		// expected
	}

	pass(name);
}

void test_poly_mul_naive() {
	const std::string name = "test_poly_mul_naive";

	{
		ivector f = {1, 2}; // 1 + 2x
		ivector g = {3, 4}; // 3 + 4x
		uint64_t m = 3;
		ivector expect = {3, 10, 8, 0, 0, 0};
		ivector got = poly_mul_naive(f, g, m);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " basic", &got, &expect);
			return;
		}
	}

	{
		// example:
		// f = 1 + x^3, g = 1 + x^2
		// product = [1,0,1,1,0,1] => reduce mod x^4 + 1 => [1, -1, 1, 1]
		ivector f = {1, 0, 0, 1};
		ivector g = {1, 0, 1};
		uint64_t m = 2;
		ivector expect = {1, -1, 1, 1};
		ivector got = poly_mul_naive(f, g, m);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " wrap", &got, &expect);
			return;
		}
	}

	{
		ivector f = {5};
		ivector g = {7};
		uint64_t m = 5;
		ivector expect = ivector(2 * m, 0);
		expect[0] = 35;
		ivector got = poly_mul_naive(f, g, m);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " const", &got, &expect);
			return;
		}
	}

	pass(name);
}

void test_poly_pow() {
	const std::string name = "test_poly_pow";

	// constant base: {2}^10 = {1024}
	{
		ivector base = {2};
		uint64_t exp = 10;
		uint64_t m = 4;
		ivector expect(2 * m, 0);
		expect[0] = 1024;
		ivector got = poly_pow(base, exp, m);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " const_pow", &got, &expect);
			return;
		}
	}

	// (1 + x)^3 = 1 + 3x + 3x^2 + x^3
	{
		ivector base = {1, 1};
		uint64_t exp = 3;
		uint64_t m = 10;
		ivector expect(2 * m, 0);
		expect[0] = 1;
		expect[1] = 3;
		expect[2] = 3;
		expect[3] = 1;
		ivector got = poly_pow(base, exp, m);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " binomial", &got, &expect);
			return;
		}
	}

	pass(name);
}

void test_poly_div_scalar() {
	const std::string name = "test_poly_div_scalar";

	// exact division
	{
		ivector p = {4, 8, 12};
		int64_t scalar = 4;
		ivector expect = {1, 2, 3};
		ivector got = poly_div_scalar(p, scalar);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " exact", &got, &expect);
			return;
		}
	}

	// non-exact division should throw
	try {
		ivector p = {3, 2};
		int64_t scalar = 2;
		(void)poly_div_scalar(p, scalar);
		fail_expectation(name + " non-exact (no exception thrown)", nullptr, nullptr);
		return;
	} catch (const std::invalid_argument &) {
		// expected
	}

	// division by zero
	try {
		ivector p = {0};
		(void)poly_div_scalar(p, 0);
		fail_expectation(name + " div-zero (no exception thrown)", nullptr, nullptr);
		return;
	} catch (const std::invalid_argument &) {
		// expected
	}

	pass(name);
}

void test_monomial_from_exponent() {
	const std::string name = "test_monomial_from_exponent";

	// exp < phi_deg
	{
		int64_t exp = 2;
		uint64_t phi = 6;
		ivector expect(phi, 0);
		expect[2] = 1;
		ivector got = monomial_from_exponent(exp, phi);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " small_exp", &got, &expect);
			return;
		}
	}

	// exp >= phi_deg (wrap odd => negative)
	{
		int64_t exp = 7; // wraps = 1, reduced = 1 => coefficient -1 at 1
		uint64_t phi = 6;
		ivector expect(phi, 0);
		expect[1] = -1;
		ivector got = monomial_from_exponent(exp, phi);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " wrap_neg", &got, &expect);
			return;
		}
	}

	// exp >= phi_deg (wrap even => positive)
	{
		int64_t exp = 12; // wraps = 2, reduced = 0 => coefficient +1 at 0
		uint64_t phi = 6;
		ivector expect(phi, 0);
		expect[0] = 1;
		ivector got = monomial_from_exponent(exp, phi);
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " wrap_pos", &got, &expect);
			return;
		}
	}

	pass(name);
}

void test_propagate_carries() {
	const std::string name = "test_propagate_carries";

	// positive carry propagation
	{
		ivector h = {15, 9};
		int64_t base = 10;
		ivector got = propagate_carries(h, base);
		ivector expect = {5, 0, 1}; // 15 -> 5 carry 1; 9 + 1 -> 10 -> 0 carry 1
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " positive", &got, &expect);
			return;
		}
	}

	// negative intermediate value
	{
		ivector h = {25, -7};
		int64_t base = 10;
		ivector got = propagate_carries(h, base);
		// manual computation:
		// start: [25, -7]
		// final: [5, 5, -1]
		ivector expect = {5, 5, -1};
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " negative-borrow", &got, &expect);
			return;
		}
	}

	// no-op when all digits < base and >= 0
	{
		ivector h = {3, 4, 5};
		int64_t base = 10;
		ivector got = propagate_carries(h, base);
		ivector expect = {3, 4, 5};
		if (!eq_ivector(got, expect)) {
			fail_expectation(name + " noop", &got, &expect);
			return;
		}
	}

	pass(name);
}

void test_FFT_linearity() {
	const std::string name = "test_FFT_linearity";

	const size_t n = 8;
	const uint64_t m = 2;
	const uint64_t phi_deg = 2 * m;

	poly_ivector f(n, ivector(phi_deg, 0));
	poly_ivector g(n, ivector(phi_deg, 0));
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < phi_deg; ++j) {
			f[i][j] = (int64_t)(i + 1 + j);
			g[i][j] = (int64_t)(2 * i - (int)j);
		}

	ivector w = monomial_from_exponent(1, phi_deg);

	try {
		poly_ivector Ff = FFT(f, w, m);
		poly_ivector Fg = FFT(g, w, m);

		// compute FFT(f+g)
		poly_ivector sum_fg(n);
		for (size_t i = 0; i < n; ++i)
			sum_fg[i] = poly_add(f[i], g[i]);
		poly_ivector Fsum = FFT(sum_fg, w, m);

		// compute Ff + Fg
		poly_ivector Ff_plus_Fg(n);
		for (size_t i = 0; i < n; ++i)
			Ff_plus_Fg[i] = poly_add(Ff[i], Fg[i]);

		if (!eq_poly_ivector(Fsum, Ff_plus_Fg)) {
			fail_expectation(name + " linearity violated", nullptr, nullptr);
			return;
		}
	} catch (const std::exception &ex) {
		fail_expectation(name + " exception: " + ex.what(), nullptr, nullptr);
		return;
	}

	pass(name);
}

void test_FFT_IFFT_large_roundtrip() {
	const std::string name = "test_FFT_IFFT_large_roundtrip";

	// choose n = 256 -> m = n/4 = 64 (so 4*m == n and x is an n-th root of unity)
	const size_t n = 256;
	const uint64_t m = 64;
	const uint64_t phi_deg = 2 * m;

	poly_ivector f(n, ivector(phi_deg, 0));
	// fill with moderate integers (not tiny, not insane)
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < phi_deg; ++j)
			f[i][j] = static_cast<int64_t>((int64_t)(i + 1) * ((int64_t)j - 7) + (int64_t)j);

	ivector w = monomial_from_exponent(1, phi_deg);

	try {
		poly_ivector F = FFT(f, w, m);
		poly_ivector H = IFFT(F, w, m);

		if (!eq_poly_ivector(H, f)) {
			ivector H_size(H.size());
			ivector f_size(f.size());
			fail_expectation(name + " large roundtrip mismatch", &H_size, &f_size);
			return;
		}
	} catch (const std::exception &ex) {
		fail_expectation(name + " exception: " + ex.what(), nullptr, nullptr);
		return;
	}

	pass(name);
}

void test_FFT_IFFT_exceptions() {
	const std::string name = "test_FFT_IFFT_exceptions";
	bool ok = true;

	// n = 3 (not a power of two)
	poly_ivector badf(3, ivector(2, 0));
	ivector w = monomial_from_exponent(1, 2);

	try {
		(void)FFT(badf, w, 1);
		fail_expectation(name + " did not throw for n=3: ", nullptr, nullptr);
		ok = false;
	} catch (const std::invalid_argument &) {
		// expected
	} catch (const std::exception &ex) {
		fail_expectation(name + " exception: " + ex.what(), nullptr, nullptr);
		ok = false;
	}

	try {
		(void)IFFT(badf, w, 1);
		fail_expectation(name + " did not throw for n=3: ", nullptr, nullptr);
		ok = false;
	} catch (const std::invalid_argument &) {
		// expected
	} catch (const std::exception &ex) {
		fail_expectation(name + " exception: " + ex.what(), nullptr, nullptr);
		ok = false;
	}

	if (ok)
		pass(name);
}

// NOTE: Only as a utility, not a real test:
static poly_ivector naive_PWC_poly(const poly_ivector &f, const poly_ivector &g, uint64_t m) {
	if (f.size() != g.size())
		throw std::invalid_argument("naive_PWC_poly: size mismatch");
	size_t n = f.size();
	uint64_t phi_deg = 2 * m;
	poly_ivector res(n, ivector(phi_deg, 0));

	for (size_t l = 0; l < n; ++l) {
		ivector accum(phi_deg, 0);
		for (size_t i = 0; i < n; ++i) {
			size_t j = (l + n - i) % n; // j such that i + j ≡ l (mod n)
			ivector prod = poly_mul_naive(f[i], g[j], m);
			accum = poly_add(accum, prod);
		}
		res[l] = accum;
	}
	return res;
}

// NOTE: Only as a utility, not a real test:
static poly_ivector naive_NWC_poly(const poly_ivector &f, const poly_ivector &g, uint64_t m) {
	if (f.size() != g.size())
		throw std::invalid_argument("naive_NWC_poly: size mismatch");
	size_t n = f.size();
	uint64_t phi_deg = 2 * m;

	// full_conv[t] is a polynomial in R (ivector length phi_deg)
	std::vector<ivector> full_conv(2 * n - 1, ivector(phi_deg, 0));

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			size_t t = i + j;
			ivector prod = poly_mul_naive(f[i], g[j], m);
			full_conv[t] = poly_add(full_conv[t], prod);
		}
	}

	poly_ivector result(n, ivector(phi_deg, 0));
	// first n contributions
	for (size_t i = 0; i < n; ++i)
		result[i] = full_conv[i];
	// subtract the tail
	for (size_t t = n; t < 2 * n - 1; ++t) {
		size_t target = t - n;
		result[target] = poly_sub(result[target], full_conv[t]);
	}
	return result;
}

void test_fast_NWC_vs_naive_poly() {
	const std::string name = "test_fast_NWC_vs_naive_poly";

	std::mt19937_64 rng(0xC0FFEE1234ULL);
	std::uniform_int_distribution<int64_t> dist(-100, 100);

	std::vector<size_t> sizes = {4, 8, 16};
	for (size_t n : sizes) {
		for (int trial = 0; trial < 20; ++trial) {
			ivector a(n), b(n);
			for (size_t i = 0; i < n; ++i) {
				a[i] = dist(rng);
				b[i] = dist(rng);
			}

			try {
				ivector naive = naive_NWC(a, b);
				ivector fast = fast_NWC(a, b);
				if (!eq_ivector(naive, fast)) {
					fail_expectation(std::format("{} did mismatch for n={}, trial={}: ", name, n, trial), &fast, &naive);
					return;
				}
			} catch (const std::exception &ex) {
				fail_expectation(std::format("{} exception for n={}, trial={}: {}", name, n, trial, ex.what()), nullptr, nullptr);
				return;
			}
		}
	}

	pass(name);
}

void test_PWC_with_PROU_vs_naive_poly() {
	const std::string name = "test_PWC_with_PROU_vs_naive_poly";

	// choose n = 8, m = 2 (phi_deg = 4) so that x is an n-PROU (4*m == n)
	const size_t n = 8;
	const uint64_t m = 2;
	const uint64_t phi_deg = 2 * m;

	poly_ivector f(n, ivector(phi_deg, 0)), g(n, ivector(phi_deg, 0));
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < phi_deg; ++j) {
			f[i][j] = static_cast<int64_t>((int64_t)(i + 1) * (int64_t)(j + 2));
			g[i][j] = static_cast<int64_t>((int64_t)(2 * i + 3) - (int64_t)j);
		}

	ivector w = monomial_from_exponent(1, phi_deg); // x is n-PROU here

	try {
		poly_ivector expected = naive_PWC_poly(f, g, m);
		poly_ivector got = PWC_with_PROU(f, g, w, m);

		if (!eq_poly_ivector(expected, got)) {
			fail_expectation(std::format("{} did mismatch for n={}, trial={}: ", name, n, m), &got[0], &expected[0]);
			return;
		}
	} catch (const std::exception &ex) {
		fail_expectation(name + " exception: " + ex.what(), nullptr, nullptr);
		return;
	}

	pass(name);
}

void test_NWC_with_PROU_vs_naive_poly() {
	const std::string name = "test_NWC_with_PROU_vs_naive_poly";

	const size_t n = 8;
	const uint64_t m = n / 2; // m = 4, phi_deg = 8; x is 2n-PROU
	const uint64_t phi_deg = 2 * m;

	poly_ivector f(n, ivector(phi_deg, 0)), g(n, ivector(phi_deg, 0));
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < phi_deg; ++j) {
			// Use small numbers but ensure some non-zero higher coefficients to exercise folding
			f[i][j] = static_cast<int64_t>((int64_t)(i + 2) * ((int64_t)j - 3));
			g[i][j] = static_cast<int64_t>((int64_t)(i + 1) + (int64_t)j % 5);
		}

	// w must be a 2n-PROU; with m = n/2, x is a 2n-PROU, so choose monomial degree 1.
	ivector w = monomial_from_exponent(1, phi_deg);

	try {
		poly_ivector expected = naive_NWC_poly(f, g, m);
		poly_ivector got = NWC_with_PROU(f, g, w, m);

		if (!eq_poly_ivector(expected, got)) {
			fail_expectation(std::format("{} did mismatch for n={}, trial={}: ", name, n, m), &got[0], &expected[0]);
			return;
		}
	} catch (const std::exception &ex) {
		fail_expectation(name + " exception: " + ex.what(), nullptr, nullptr);
		return;
	}

	pass(name);
}

void test_naive_NWC_basic_example() {
	const std::string name = "test_naive_NWC_basic_example";

	ivector f = {1, 2, 3};
	ivector g = {4, 5, 6};
	ivector expect = {-23, -5, 28}; // computed by hand/full_conv folding
	ivector got = naive_NWC(f, g);

	if (!eq_ivector(got, expect)) {
		fail_expectation(name, &got, &expect);
		return;
	}

	pass(name);
}

void test_txt_input_small() {
	const std::string name = "test_txt_input_small";

	std::string a_str, b_str, ab_str;
	read_txt(a_str, "../large_numbers.txt", 0);
	read_txt(b_str, "../large_numbers.txt", 1);
	read_txt(ab_str, "../large_numbers.txt", 2);

	ivector a_vector, b_vector, ab_vector_correct;
	string_to_ivector(a_str, a_vector);
	string_to_ivector(b_str, b_vector);
	string_to_ivector(ab_str, ab_vector_correct);

	ivector ab_vector_given = propagate_carries(full_product(a_vector, b_vector), (1ULL << BITS_PER_PART));

	if (!eq_ivector(ab_vector_correct, ab_vector_given)) {
		fail_expectation(name + " did mismatch", &ab_vector_correct, &ab_vector_given);
		return;
	}

	pass(name);
}

void test_txt_input_large() {
	const std::string name = "test_txt_input_large";

	std::string a_str, b_str, ab_str;
	read_txt(a_str, "../large_numbers.txt", 4);
	read_txt(b_str, "../large_numbers.txt", 5);
	read_txt(ab_str, "../large_numbers.txt", 6);

	ivector a_vector, b_vector, ab_vector_correct;
	string_to_ivector(a_str, a_vector);
	string_to_ivector(b_str, b_vector);
	string_to_ivector(ab_str, ab_vector_correct);

	ivector ab_vector_given = propagate_carries(full_product(a_vector, b_vector), (1ULL << BITS_PER_PART));

	std::cout << "Correct product:\n";
	for (size_t i = 0; i < 10 && i < ab_vector_correct.size(); i++) {
		std::cout << static_cast<int64_t>(ab_vector_correct[i]) << std::endl;
	}
	std::cout << "Given product:\n";
	for (size_t i = 0; i < 10 && i < ab_vector_given.size(); i++) {
		std::cout << static_cast<int64_t>(ab_vector_given[i]) << std::endl;
	}

	if (!eq_ivector(ab_vector_correct, ab_vector_given)) {
		fail_expectation(name + " did mismatch", &ab_vector_correct, &ab_vector_given);
		return;
	}

	pass(name);
}
