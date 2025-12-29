#include "convolution.hpp"
#include "utils.hpp"

using namespace std;

int main() {
	string a_str, b_str, ab_str;
	read_txt(a_str, "../large_numbers.txt", 0);
	read_txt(b_str, "../large_numbers.txt", 1);
	read_txt(ab_str, "../large_numbers.txt", 2);

	ivector a_vector, b_vector, ab_vector_correct;
	string_to_ivector(a_str, a_vector);
	string_to_ivector(b_str, b_vector);
	string_to_ivector(ab_str, ab_vector_correct);

	cout << "Given `a`:\n";
	print_vector(a_vector);

	cout << "Given `b`:\n";
	print_vector(b_vector);

	cout << a_vector.size() << " " << b_vector.size() << endl;

	ivector ab_vector_given = fast_NWC(a_vector, b_vector);

	cout << "Correct product:\n";
	print_vector(ab_vector_correct);
	cout << "Given product:\n";
	print_vector(ab_vector_given);

	return 0;
}
