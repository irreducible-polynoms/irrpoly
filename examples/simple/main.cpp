// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include <irrpoly.h>

#include <iostream>

using namespace irrpoly;

auto main() -> int {
    std::cout << std::noboolalpha;

    auto field = make_gf(3);
    auto poly = gfpoly(field, {2, 1, 1});

    std::cout << poly << " "
              << is_irreducible(poly) << " "
              << is_primitive(poly) << std::endl;

    return 0;
}
