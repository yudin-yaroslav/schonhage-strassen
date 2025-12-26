#include <complex>
#include <iostream>
#include <numbers>
#include <cmath>
#include <vector>

#define DEBUG_MULT false

using namespace std;
using namespace std::complex_literals;
using namespace std::numbers;

typedef complex<double> dcomplex;
typedef vector<double> dvector;

// Helper functions
void print_polynomial(const dvector &coeffs) {
	for (int i = coeffs.size() - 1; i >= 0; i--) {
		if (coeffs[i] < 0) cout << "- " << -coeffs[i] << "x^" << i << " ";
		else if (coeffs[i] > 0) cout << "+ " << coeffs[i] << "x^" << i << " ";
	}
	cout << endl;
}

template <typename T>
void print_vector(const vector<T>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        cout << vec[i] << endl;
    }
}

bool approx_equal(dcomplex a, dcomplex b, double eps = 1e-9) {
    return abs(a.real() - b.real()) < eps
        && abs(a.imag() - b.imag()) < eps;
}


// FFT and IFFT functions
vector<dcomplex> FFT(dvector coeffs, int n) {
	if ((n & (n - 1)) != 0) {
		throw "FFT: The length of the coefficient list is not a power of 2";
	}
	if (coeffs.size() > n) {
		throw "FFT: Some coefficients were truncated";
	}
	if (n == 1) {
		return {dcomplex(coeffs[0], 0.0)};
	}

	double w_angle = 2 * pi / n;
	dcomplex w_n(cos(w_angle), sin(w_angle));

	coeffs.resize(n);

	dvector coeffs_even(n/2);
	dvector coeffs_odd(n/2);

	for (int i=0; i<n/2; i++) {
		coeffs_even[i] = coeffs[2*i];
		coeffs_odd[i] = coeffs[2*i + 1];
	}

	vector<dcomplex> values_even = FFT(coeffs_even, n/2);
	vector<dcomplex> values_odd = FFT(coeffs_odd, n/2);

	vector<dcomplex> values(n);
	
	dcomplex w(1.0, 0.0);
	for (int i=0; i<n/2; i++) {
		values[i] = values_even[i] + w*values_odd[i];
		values[i + n/2] = values_even[i] - w*values_odd[i];
		w *= w_n;
	}

	return values;
}

vector<dcomplex> IFFT_complex(vector<dcomplex> values, int n) {
	if ((n & (n - 1)) != 0) {
		throw "IFFT: The length of the coefficient list is not a power of 2";
	}
	if (values.size() > n) {
		throw "IFFT: Some values were truncated";
	}
	if (n == 1) {
		return {values[0]};
	}

	double w_angle = -2 * pi / n;
	dcomplex w_n(cos(w_angle), sin(w_angle));

	values.resize(n);

	vector<dcomplex> values_even(n/2);
	vector<dcomplex> values_odd(n/2);

	for (int i=0; i<n/2; i++) {
		values_even[i] = values[2*i];
		values_odd[i] = values[2*i + 1];
	}

	vector<dcomplex> coeffs_even = IFFT_complex(values_even, n/2);
	vector<dcomplex> coeffs_odd = IFFT_complex(values_odd, n/2);

	vector<dcomplex> coeffs(n);
	
	dcomplex w(1.0, 0.0);
	for (int i=0; i<n/2; i++) {
		coeffs[i] = (coeffs_even[i] + w*coeffs_odd[i]) / 2.0;
		coeffs[i + n/2] = (coeffs_even[i] - w*coeffs_odd[i]) / 2.0;
		w *= w_n;
	}

	return coeffs;
}

dvector IFFT(vector<dcomplex> values, int n) {
	vector<dcomplex> result_complex = IFFT_complex(values, n);

	dvector result_real(n);
	for (int i=0; i<n; i++) {
		result_real[i] = result_complex[i].real();
	}

	return result_real;
}

// Fast polynomial multiplication function
dvector mult_polynomial(dvector A_coeffs, dvector B_coeffs) {
	int A_degree = A_coeffs.size() - 1;
	int B_degree = B_coeffs.size() - 1;

	int n = 1;
	while (n < A_degree + B_degree + 1) n <<= 1;

	vector<dcomplex> A_values = FFT(A_coeffs, n);
	vector<dcomplex> B_values = FFT(B_coeffs, n);

	if (DEBUG_MULT) {
		cout << "A_values:" << endl;
		print_vector(A_values);
		cout << "B_values:" << endl;
		print_vector(B_values);
	}

	vector<dcomplex> C_values(n);
	for (int i=0; i<n; i++) {
		C_values[i] = A_values[i] * B_values[i];
	}
	if (DEBUG_MULT) {
		cout << "C_values:" << endl;
		print_vector(C_values);
	}

	dvector C_coeffs = IFFT(C_values, n);
	if (DEBUG_MULT) {
		cout << "C_coeffs:" << endl;
		print_polynomial(C_coeffs);
	}

	return C_coeffs;
}

// Testing functions
bool test_FFT() {
	dvector coeffs = { 3, -2, 5, 1, -4 };

	vector<dcomplex> values_correct = {
	    { 3.0000000000,  0.0000000000 },
	    { 4.8786796564, 4.2928932180 },
	    { -6.0000000000,  -3.0000000000 },
	    { 9.1213203436,  -5.7071067812 },
	    { 5.0000000000,  0.0000000000 },
	    { 9.1213203436, 5.7071067812 },
	    { -6.0000000000, 3.0000000000 },
	    { 4.8786796564,  -4.2928932180 }
	};
	vector<dcomplex> values_answer = FFT(coeffs, 8);

    for (int i = 0; i < values_correct.size(); i++) {
		if (!approx_equal(values_answer[i], values_correct[i])) {
			return false;
        }
    }
	return true;
}

bool test_IFFT() {
    dvector coeffs = { -1, 0, 2, -3, 4 };
    int n = 8;

	coeffs.resize(n, 0.0);

    vector<dcomplex> values = FFT(coeffs, n);
    dvector coeffs_new = IFFT(values, n);

    for (int i = 0; i < n; i++) {
        if (!approx_equal(coeffs_new[i], coeffs[i])) {
            return false;
        }
    }
    return true;
}

bool test_mult_polynomial() {
	dvector A_coeffs = { 3, -2, 5, 1, -4 };
	dvector B_coeffs = { -1, 4, 2, -3 };

	dvector C_coeffs_correct = { -3, 14, -7, 6, 24, -29, -11, 12 };
	dvector C_coeffs_answer = mult_polynomial(A_coeffs, B_coeffs);

	if (DEBUG_MULT) {
		print_vector(C_coeffs_answer);
	}

    for (int i = 0; i < 8; i++) {
        if (!approx_equal(C_coeffs_answer[i], C_coeffs_correct[i])) {
            return false;
        }
    }
	return true;
}


int main() {
	dvector A_coeffs;
	dvector B_coeffs;
	dvector C_coeffs;

	bool test_FFT_result = test_FFT();
	bool test_IFFT_result = test_IFFT();
	bool test_mult_polynomial_result = test_mult_polynomial();

	cout << "Test FFT: " << test_FFT_result << endl;
	cout << "Test IFFT (inversibility): " << test_IFFT_result << endl;
	cout << "Test fast polynomial multiplication: " << test_mult_polynomial_result << endl;
}
