#include "utils.hpp"
#include <algorithm>
#include <fstream>

// Data loader
void string_to_ivector(std::string str, ivector &vec) {
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

void read_txt(std::string &out, const std::string &filename, size_t line_index) {
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

// Polynomial arithmetics
void poly_add(ivector &dst, const ivector &f, const ivector &g) {
	if (f.size() != g.size())
		throw std::invalid_argument("poly_add: Given `f` and `g` have different sizes");
	dst.resize(f.size());
	std::transform(f.begin(), f.end(), g.begin(), dst.begin(), std::plus<int64_t>{});
}

void poly_sub(ivector &dst, const ivector &f, const ivector &g) {
	if (f.size() != g.size())
		throw std::invalid_argument("poly_sub: Given `f` and `g` have different sizes");
	dst.resize(f.size());
	std::transform(f.begin(), f.end(), g.begin(), dst.begin(), std::minus<int64_t>{});
}

// WARNING: R = Z[x]/(x^{2m} + 1)
void poly_mul(ivector &dst, const ivector &f, const ivector &g, uint64_t m) {
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

	if (dst.size() > 2 * m)
		throw std::invalid_argument("poly_mul: Given `dst` doesn't have the correct size of at most 2*m");

	dst.assign(2 * m, 0);

	for (size_t i = 0; i < prod.size(); ++i) {
		size_t q = i / phi_deg;
		size_t r = i % phi_deg;
		int64_t sign = (q % 2 == 0) ? 1 : -1;
		dst[r] += sign * prod[i];
	}
}

void poly_pow(ivector &dst, const ivector &base, uint64_t exp, uint64_t m) {
	const size_t phi = 2 * m;

	ivector cur(phi, 0);
	for (size_t i = 0; i < base.size(); ++i) {
		size_t r = i % phi;
		size_t q = i / phi;
		cur[r] += (q & 1) ? -base[i] : base[i];
	}

	ivector res(phi, 0);
	res[0] = 1;
	ivector tmp;

	while (exp) {
		if (exp & 1) {
			poly_mul(tmp, res, cur, m);
			res.swap(tmp);
		}
		exp >>= 1;
		if (exp) {
			poly_mul(tmp, cur, cur, m);
			cur.swap(tmp);
		}
	}

	dst.swap(res);
}

void poly_div_scalar(ivector &f, int64_t scalar) {
	if (scalar == 0)
		throw std::invalid_argument("poly_div_scalar: Division by zero in poly_div_scalar");

	for (size_t i = 0; i < f.size(); ++i) {
		if ((f[i] % scalar) != 0) {
			throw std::invalid_argument("poly_div_scalar: Non-exact division when computing IFFT (coeff not divisible by n)");
		}
		f[i] = f[i] / scalar;
	}
}

ivector poly_mul_by_x_power(const ivector &p, uint64_t power, uint64_t m) {
	const size_t phi = 2 * m;
	if (p.size() != phi) {
		throw std::invalid_argument("mul_by_x_power: Input `p` must be padded to size 2*m");
	}
	uint64_t r = power % phi;
	int overall_sign = ((power / phi) % 2 == 0) ? 1 : -1;
	ivector res(phi, 0);
	for (size_t j = 0; j < phi; ++j) {
		if (p[j] == 0)
			continue;
		size_t new_pos = j + r;
		int local_sign = overall_sign;
		if (new_pos >= phi) {
			new_pos -= phi;
			local_sign = -local_sign;
		}
		res[new_pos] += local_sign * p[j];
	}
	return res;
}
