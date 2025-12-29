#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

using ivector = std::vector<uint64_t>;
using poly_ivector = std::vector<ivector>;

constexpr uint8_t HEX_DIGITS_PER_WORD = 32;
constexpr uint8_t NIBBLE_SHIFT = 4;
constexpr uint8_t BITS_PER_PART = HEX_DIGITS_PER_WORD * NIBBLE_SHIFT;

// Debug
template <typename T> void print_polynomial(const std::vector<T> &coeffs) {
	for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i) {
		if (coeffs[i] < 0)
			std::cout << "- " << -coeffs[i] << "x^" << i << " ";
		else if (coeffs[i] > 0)
			std::cout << "+ " << coeffs[i] << "x^" << i << " ";
	}
	std::cout << std::endl;
}

template <typename T> void print_vector(const std::vector<T> &vec) {
	for (const auto &x : vec)
		std::cout << x << std::endl;
}

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
		vec.push_back(value);
	}

	std::reverse(vec.begin(), vec.end());
}

inline void read_txt(std::string &out, const std::string &filename, size_t line_index) {
	std::ifstream file(filename);
	file.exceptions(std::ifstream::badbit);

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

// Modular arithmetics (slow)
inline ivector poly_add(const ivector &f, const ivector &g) {
	if (f.size() != g.size()) {
		throw std::invalid_argument("poly_add: Given `f` and `g` have different sizes");
	}
	const uint64_t n = f.size();

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
	const uint64_t n = f.size();

	ivector res(n);
	for (size_t i = 0; i < n; ++i) {
		res[i] = f[i] - g[i];
	}

	return res;
}

inline ivector poly_mul_naive(const ivector &f, const ivector &g, uint64_t m) {
	// WARNING: R = Z[x]/(x^{2m} + 1)
	uint64_t phi_deg = 2 * m;

	size_t deg_f = f.size() ? f.size() - 1 : 0;
	size_t deg_g = g.size() ? g.size() - 1 : 0;
	size_t deg_prod = deg_f + deg_g;

	std::vector<__int128> prod(deg_prod + 1, 0);
	for (size_t i = 0; i < deg_f; ++i) {
		for (size_t j = 0; j < deg_g; ++j) {
			prod[i + j] += static_cast<__int128>(f[i]) * static_cast<__int128>(g[j]);
		}
	}

	ivector res(phi_deg, 0);
	for (size_t i = 0; i < prod.size(); ++i) {
		size_t q = i / phi_deg;
		size_t r = i % phi_deg;
		int64_t sign = (q % 2 == 0) ? 1 : -1;
		res[r] += static_cast<int64_t>(sign * prod[i]);
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

inline ivector poly_div_scalar(const ivector &p, uint64_t scalar) {
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
