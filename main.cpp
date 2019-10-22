#include <iostream>

#include "polynomialgf.hpp"

using namespace std;

int main() {
    // irreducible
    polynomialgf<2> p0({ 1, 1, 1 });
    cout << "expected: 1, got: " << is_irreducible(p0) << endl;

    // primitive
    polynomialgf<2> p1({ 1, 0, 1, 0, 0, 1 });
    cout << "expected: 1, got: " << is_primitive(p1) << endl;
    polynomialgf<3> p2({ 2, 1, 0, 2, 1, 0, 0, 0, 1 });
    cout << "expected: 1, got: " << is_primitive(p2) << endl;
    polynomialgf<5> p3({ 2, 2, 1, 0, 1 });
    cout << "expected: 1, got: " << is_primitive(p3) << endl;

    //non-primitive
    polynomialgf<3> p4({ 1, 1, 2 });
    cout << "expected: 0, got: " << is_primitive(p4) << endl;

    return 0;
}
