#include <irrpoly.h>

#include <iostream>

using namespace irrpoly;

auto main() -> int {
    std::cout << "TEST DATA" << std::endl;
    std::cout << std::noboolalpha;
    for (auto base : {2, 3, 5, 7}) {
        auto gf = make_gf(base);
        for (size_t i = 4; i < 10; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                gfpoly poly = gfpoly::random(gf, i);
                std::cout << is_irreducible_benor(poly) << " "
                          << is_primitive_definition(poly) << " "
                          << base << " " << poly << std::endl;
            }
        }
    }
    return 0;
}
