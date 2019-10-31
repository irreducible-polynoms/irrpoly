#include <iostream>

#include "polynomialgf.hpp"

using namespace std;

int main() {
//    // irreducible
//    polynomialgf<2> p0({1, 1, 1});
//    cout << "expected: 1, got: " << is_irreducible_berlekamp(p0) << endl;
//    cout << "expected: 1, got: " << is_irreducible_rabin(p0) << endl;
//
//    // primitive
//    polynomialgf<2> p1({1, 0, 1, 0, 0, 1});
//    cout << "expected: 1, got: " << is_primitive(p1) << endl;
//    polynomialgf<3> p2({2, 1, 0, 2, 1, 0, 0, 0, 1});
//    cout << "expected: 1, got: " << is_primitive(p2) << endl;
//    polynomialgf<5> p3({2, 2, 1, 0, 1});
//    cout << "expected: 1, got: " << is_primitive(p3) << endl;
//
//    //non-primitive
//    polynomialgf<2> p5({1, 1});
//    cout << "expected: 0, got: " << is_primitive(p5) << endl;
//    polynomialgf<2> p6({1, 1, 1, 1, 1});
//    cout << "expected: 0, got: " << is_primitive(p6) << endl;
//    polynomialgf<2> p7({1, 0, 0, 1, 0, 0, 1});
//    cout << "expected: 0, got: " << is_primitive(p7) << endl;
//    polynomialgf<2> p8({1, 1, 0, 0, 1, 0, 1});
//    cout << "expected: 0, got: " << is_primitive(p8) << endl;
//    polynomialgf<2> p9({1, 1, 1, 0, 1, 0, 1});
//    cout << "expected: 0, got: " << is_primitive(p9) << endl;

    // TESTING
    const uint32_t P = 5; // change P to what ever you want
    for (size_t i = 4; i < 24; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            polynomialgf<P> p = random<P>(i);
            cout << (is_irreducible_rabin(p) ? '1' : '0') << " " // also is_irreducible_berlekamp
                 << (is_primitive(p) ? '1' : '0') << " " << P << " "
                 << p << endl;
        }
    }
    // copy output and save to file named "input" (WITHOUT EXTENSION!)
    // use Wolfram Mathematica
    // define function for validating the result
    // f[irr_, prim_, m_, p_] :=
    // irr == Boole[
    //    IrreduciblePolynomialQ[Dot[Power[x, Range[0, Length[p] - 1]], p],
    //     Modulus -> m]] &&
    //  prim == Boole[
    //    PrimitivePolynomialQ[Dot[Power[x, Range[0, Length[p] - 1]], p],
    //     Modulus -> m]]
    // open "input" file for read (CHANGE PATH TO FILE!)
    // file = OpenRead["/Users/vadimpiven/Downloads/input"]
    // validate results line by line
    // While[Not[(irr = Read[file, Number]) === EndOfFile],
    // prim = Read[file, Number]; m = Read[file, Number];
    // p = Read[StringToStream[#], Number] & /@
    //   StringSplit[StringTrim[ReadLine[file], ("{" | "}" | " ") ...],
    //    ", "]; Print[f[irr, prim, m, p]]]
    // EXPECTED RESULT IS COLUMN OF "True", if there is some False - this program has some errors
    // close the file
    // Close[file]

    return 0;
}
