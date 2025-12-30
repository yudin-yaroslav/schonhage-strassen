#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

using ivector = std::vector<int64_t>;
using poly_ivector = std::vector<ivector>;

constexpr uint8_t HEX_DIGITS_PER_WORD = 4;
constexpr uint8_t NIBBLE_SHIFT = 4;
constexpr uint8_t BITS_PER_PART = HEX_DIGITS_PER_WORD * NIBBLE_SHIFT;

// Helper utils
static bool eq_ivector(const ivector &a, const ivector &b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); ++i)
		if (a[i] != b[i])
			return false;
	return true;
}

static bool eq_poly_ivector(const poly_ivector &a, const poly_ivector &b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); ++i)
		if (!eq_ivector(a[i], b[i]))
			return false;
	return true;
}

static void print_ivector(const ivector &v) {
	std::cout << "{";
	for (size_t i = 0; i < v.size(); ++i) {
		if (i)
			std::cout << ", ";
		std::cout << v[i];
	}
	std::cout << "}";
}

static void fail_expectation(const std::string &testname, ivector *got, ivector *expect) {
	std::cout << "\033[31m" << "FAIL " << testname << "\033[0m\t\n";
	if (got != nullptr && expect != nullptr) {
		std::cout << "Got:\n";
		print_ivector(*got);
		std::cout << "\nExpected:\n";
		print_ivector(*expect);
		std::cout << std::endl;
	}
}

static void pass(const std::string &testname) { std::cout << "\033[32m" << "PASS " << testname << "\033[0m\t\n"; }

// Data loader
inline uint32_t hexval(char c) {
	if ('0' <= c && c <= '9')
		return c - '0';
	if ('a' <= c && c <= 'f')
		return c - 'a' + 10;
	if ('A' <= c && c <= 'F')
		return c - 'A' + 10;
	throw std::invalid_argument("Invalid hex digit");
}

inline void string_to_ivector(std::string str, ivector &vec) {
	if (str.size() % HEX_DIGITS_PER_WORD != 0) {
		str.insert(0, HEX_DIGITS_PER_WORD - (str.size() % HEX_DIGITS_PER_WORD), '0');
	}

	vec.clear();
	vec.reserve(str.size() / HEX_DIGITS_PER_WORD);

	for (size_t i = 0; i < str.size(); i += HEX_DIGITS_PER_WORD) {
		uint64_t value = 0;
		for (size_t j = 0; j < HEX_DIGITS_PER_WORD; ++j) {
			value = (value << NIBBLE_SHIFT) | hexval(str[i + j]);
		}
		vec.push_back(static_cast<int64_t>(value));
	}

	std::reverse(vec.begin(), vec.end());
}

inline void read_txt(std::string &out, const std::string &filename, size_t line_index) {
	std::ifstream file(filename);
	file.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	try {
		std::string line;
		size_t counter = 0;
		while (getline(file, line)) {
			if (counter == line_index) {
				out = line;
				return;
			}
			counter++;
		}

		file.close();
	} catch (const std::exception &ex) {
		throw std::runtime_error(std::string("I/O error while reading file: ") + ex.what() + "\n");
	}
}

// Polynomial arithmetic
inline ivector poly_add(const ivector &f, const ivector &g) {
	if (f.size() != g.size()) {
		throw std::invalid_argument("poly_add: Given `f` and `g` have different sizes");
	}
	const size_t n = f.size();

	ivector res(n);
	for (size_t i = 0; i < n; ++i) {
		res[i] = f[i] + g[i];
	}

	return res;
}

inline ivector poly_sub(const ivector &f, const ivector &g) {
	if (f.size() != g.size()) {
		throw std::invalid_argument("poly_sub: Given `f` and `g` have different sizes");
	}
	const size_t n = f.size();

	ivector res(n);
	for (size_t i = 0; i < n; ++i) {
		res[i] = f[i] - g[i];
	}

	return res;
}

inline ivector poly_mul_naive(const ivector &f, const ivector &g, uint64_t m) {
	// WARNING: R = Z[x]/(x^{2m} + 1)
	uint64_t phi_deg = 2 * m;

	size_t deg_f = f.size() - 1;
	size_t deg_g = g.size() - 1;
	size_t deg_prod = deg_f + deg_g;

	ivector prod(deg_prod + 1, 0);
	for (size_t i = 0; i < f.size(); ++i) {
		for (size_t j = 0; j < g.size(); ++j) {
			prod[i + j] += f[i] * g[j];
		}
	}

	ivector res(phi_deg, 0);
	for (size_t i = 0; i < prod.size(); ++i) {
		size_t q = i / phi_deg;
		size_t r = i % phi_deg;
		int64_t sign = (q % 2 == 0) ? 1 : -1;
		res[r] += sign * prod[i];
	}

	return res;
}

inline ivector poly_pow(ivector base, uint64_t exp, uint64_t m) {
	ivector result = {1};
	while (exp > 0) {
		if (exp & 1)
			result = poly_mul_naive(result, base, m);
		base = poly_mul_naive(base, base, m);
		exp >>= 1;
	}
	return result;
}

inline ivector poly_div_scalar(const ivector &p, int64_t scalar) {
	if (scalar == 0)
		throw std::invalid_argument("poly_div_scalar: Division by zero in poly_div_scalar");
	ivector res = p;
	for (size_t i = 0; i < res.size(); ++i) {
		if ((res[i] % scalar) != 0) {
			throw std::invalid_argument("poly_div_scalar: Non-exact division when computing IFFT (coeff not divisible by n)");
		}
		res[i] = res[i] / scalar;
	}
	return res;
}

inline ivector monomial_from_exponent(int64_t exp, uint64_t phi_deg) {
	ivector poly(phi_deg, 0);
	if (exp < phi_deg) {
		poly[exp] = 1;
	} else {
		// Reduce modulo x^{phi_deg} + 1
		size_t wraps = exp / phi_deg;
		size_t reduced_exp = exp % phi_deg;
		poly[reduced_exp] = (wraps % 2 == 0) ? 1 : -1;
	}
	return poly;
}

inline ivector propagate_carries(ivector h, int64_t base) {
	h.push_back(0);

	for (size_t i = 0; i + 1 < h.size(); ++i) {
		int64_t carry = h[i] / base;
		h[i] %= base;
		if (h[i] < 0) {
			int64_t borrow = (std::abs(h[i]) + base - 1) / base;
			h[i] += borrow * base;
			carry -= borrow;
		}
		h[i + 1] += carry;
	}

	while (!h.empty() && h.back() >= base) {
		int64_t carry = h.back() / base;
		h.back() %= base;
		h.push_back(carry);
	}

	while (h.back() == 0) {
		h.pop_back();
	}

	return h;
}
