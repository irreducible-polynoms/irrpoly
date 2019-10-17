#include <iostream>

#include "polynomialgf.hpp"

using namespace std;

int main() {
    polynomialgf<2> test({ 1, 1, 1 }); // 1 + x + x^2 (mod 2)

    cout << is_irreducible(test) << endl; // Berlekamp test

    return 0;
}
