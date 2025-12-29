#pragma once

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

using ivector = std::vector<uint64_t>;

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
constexpr uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m) { return static_cast<uint64_t>((unsigned __int128)a * b % m); }

constexpr uint64_t add_mod(uint64_t a, uint64_t b, uint64_t m) {
	uint64_t s = a + b;
	if (s >= m || s < a)
		s -= m;
	return s;
}

constexpr uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t m) { return (a >= b) ? (a - b) : (a + m - b); }

constexpr uint64_t pow_mod(uint64_t base, uint64_t exp, uint64_t m) {
	uint64_t res = 1 % m;
	base %= m;
	while (exp) {
		if (exp & 1)
			res = mul_mod(res, base, m);
		base = mul_mod(base, base, m);
		exp >>= 1;
	}
	return res;
}

//
// // NTT prime generation
// constexpr bool is_prime(uint64_t n) {
// 	if (n < 2)
// 		return false;
// 	if ((n & 1) == 0)
// 		return n == 2;
//
// 	for (uint64_t i = 3; i * i <= n; i += 2) {
// 		if (n % i == 0)
// 			return false;
// 	}
// 	return true;
// }
//
// consteval uint64_t get_ntt_prime(uint8_t k) {
// 	const uint64_t pow2 = (1ULL << k);
// 	uint64_t n = pow2 + 1;
//
// 	while (!is_prime(n)) {
// 		n += pow2;
// 	}
//
// 	return n;
// }
//
// consteval std::array<uint64_t, 32> build_ntt_prime_table() {
// 	std::array<uint64_t, 32> table{};
// 	for (size_t i = 0; i < 32; ++i) {
// 		table[i] = get_ntt_prime(i);
// 	}
// 	return table;
// }
// inline constinit const std::array<uint64_t, 32> ntt_primes = build_ntt_prime_table();
//
// // NTT elements of order k generation
// consteval bool is_primitive_root(uint64_t g, uint64_t mod) {
// 	if (mod < 2 || g % mod == 0 || !is_prime(mod))
// 		return false;
//
// 	uint64_t phi = mod - 1;
// 	uint64_t m = phi;
//
// 	if ((m & 1) == 0) {
// 		if (pow_mod(g, phi / 2, mod) == 1)
// 			return false;
// 		while ((m & 1) == 0)
// 			m >>= 1;
// 	}
//
// 	for (uint64_t prime = 3; prime * prime <= m; prime += 2) {
// 		if (m % prime == 0) {
// 			if (pow_mod(g, phi / prime, mod) == 1)
// 				return false;
// 			while (m % prime == 0)
// 				m /= prime;
// 		}
// 	}
//
// 	if (m > 1) {
// 		if (pow_mod(g, phi / m, mod) == 1)
// 			return false;
// 	}
//
// 	return true;
// }
//
// consteval uint64_t get_ntt_order_k(uint8_t k) {
// 	const uint64_t pow2 = (1ULL << k);
// 	const uint64_t prime = ntt_primes[k];
//
// 	for (uint64_t g = 2; g < prime; ++g) {
// 		if (is_primitive_root(g, prime)) {
// 			return pow_mod(g, (prime - 1) / pow2, prime);
// 		}
// 	}
//
// 	return 0;
// }
//
// consteval std::array<uint64_t, 32> build_ntt_order_k() {
// 	std::array<uint64_t, 32> table{};
// 	for (size_t i = 0; i < 32; ++i) {
// 		table[i] = get_ntt_order_k(static_cast<uint8_t>(i));
// 	}
// 	return table;
// }
// inline constinit const std::array<uint64_t, 32> ntt_order_k = build_ntt_order_k();

inline bool is_prou(uint64_t w, uint64_t n, uint64_t mod) {
	if (pow_mod(w, n, mod) != 1)
		return false;

	for (uint64_t i = 1; i < n; ++i) {
		uint64_t sum = 1;
		uint64_t cur = 1;
		for (uint64_t j = 1; j < n; ++j) {
			cur = mul_mod(cur, pow_mod(w, i, mod), mod);
			sum = add_mod(sum, cur, mod);
		}
		if (sum != 0) {
			return false;
		}
	}

	return true;
}
