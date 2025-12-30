#include "product.hpp"
#include "tests.hpp"

using namespace std;

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
	test_txt_input_large();

	return 0;
}
