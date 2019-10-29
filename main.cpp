#include <iostream>

#include "polynomialgf.hpp"

using namespace std;

int main() {
    // irreducible
    polynomialgf<2> p0({1, 1, 1});
    cout << "expected: 1, got: " << is_irreducible_berlekamp(p0) << endl;
    cout << "expected: 1, got: " << is_irreducible_rabin(p0) << endl;

    // primitive
    polynomialgf<2> p1({1, 0, 1, 0, 0, 1});
    cout << "expected: 1, got: " << is_primitive(p1) << endl;
    polynomialgf<3> p2({2, 1, 0, 2, 1, 0, 0, 0, 1});
    cout << "expected: 1, got: " << is_primitive(p2) << endl;
    polynomialgf<5> p3({2, 2, 1, 0, 1});
    cout << "expected: 1, got: " << is_primitive(p3) << endl;

    //non-primitive
    polynomialgf<2> p5({1, 1});
    cout << "expected: 0, got: " << is_primitive(p5) << endl;
    polynomialgf<2> p6({1, 1, 1, 1, 1});
    cout << "expected: 0, got: " << is_primitive(p6) << endl;
    polynomialgf<2> p7({1, 0, 0, 1, 0, 0, 1});
    cout << "expected: 0, got: " << is_primitive(p7) << endl;
    polynomialgf<2> p8({1, 1, 0, 0, 1, 0, 1});
    cout << "expected: 0, got: " << is_primitive(p8) << endl;
    polynomialgf<2> p9({1, 1, 1, 0, 1, 0, 1});
    cout << "expected: 0, got: " << is_primitive(p9) << endl;

    return 0;
}
