#include "mult.hpp"
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

	ivector ab_vector_given = multiply(a_vector, b_vector);

	cout << "Correct: " << endl;
	for (int i = 0; i < 10; i++) {
		cout << ab_vector_correct[i] << endl;
	}
	cout << "Given: " << endl;
	for (int i = 0; i < 10; i++) {
		cout << ab_vector_given[i] << endl;
	}

	return 0;
}
