#include "tests.hpp"
#include <chrono>
#include <iostream>

using namespace std;
using namespace std::chrono;

int main() {
	test_poly_add();
	test_poly_sub();
	test_poly_mul_naive();
	test_poly_pow();
	test_poly_div_scalar();
	test_monomial_from_exponent();
	test_propagate_carries();
	test_FFT_linearity();
	test_FFT_IFFT_large_roundtrip();
	test_FFT_IFFT_exceptions();
	test_fast_NWC_vs_naive_poly();
	test_PWC_with_PROU_vs_naive_poly();
	test_NWC_with_PROU_vs_naive_poly();
	test_naive_NWC_basic_example();

	test_txt_input_small();
	const auto start = steady_clock::now();
	test_txt_input_large();
	const auto end = steady_clock::now();

	const auto ns = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	cout << "duration: " << ns / 1e6 << " ms\n";

	return 0;
}
