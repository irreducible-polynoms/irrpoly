#include <iostream>

#include "polynomialgf.hpp"

using namespace std;

int main() {
    PolynomialGF<2> test({ 1, 1, 1 }); // 1 + x + x^2 (mod 2)

    cout << isIrreducible(test) << endl; // berlekamp check

	cin.get();
    return 0;
}