#include <iostream>

#include <irrpoly/gfcheck.hpp>

using namespace irrpoly;

int main() {
    std::cout << "TEST DATA" << std::endl;
    for (auto base : {2, 3, 5, 7}) {
        auto gf = make_gf(base);
        for (size_t i = 4; i < 10; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                gfpoly poly = gfpoly::random(gf, i);
                std::cout << (is_irreducible_rabin(poly) ? '1' : '0') << " "
                          << (is_primitive_definition(poly) ? '1' : '0') << " " << base << " "
                          << poly << std::endl;
            }
        }
    }
    return 0;
}
