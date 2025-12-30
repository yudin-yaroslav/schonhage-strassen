#pragma once

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

void string_to_ivector(std::string str, ivector &vec);
void read_txt(std::string &out, const std::string &filename, size_t line_index);

// Polynomial arithmetic
void poly_add(ivector &dst, const ivector &f, const ivector &g);
void poly_sub(ivector &dst, const ivector &f, const ivector &g);
void poly_mul(ivector &dst, const ivector &f, const ivector &g, uint64_t m);
void poly_pow(ivector &dst, const ivector &base, uint64_t exp, uint64_t m);
void poly_div_scalar(ivector &f, int64_t scalar);

// Other utils
inline ivector monomial_from_exponent(int64_t exp, uint64_t phi_deg) {
	ivector poly(phi_deg, 0);
	if (exp < phi_deg) {
		poly[exp] = 1;
	} else {
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

ivector poly_mul_by_x_power(const ivector &p, uint64_t power, uint64_t m);
