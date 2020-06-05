#include <irrpoly.h>

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

#include <iostream>
#include <iomanip>

using namespace irrpoly;

auto main() -> int {
    nmod_poly_t poly0;
    nmod_poly_init(poly0, 3);

    auto field = make_gf(3);
    auto poly1 = gfpoly(field, {0, 1});
    nmod_poly_set_coeff_ui(poly0, 0, 0);
    nmod_poly_set_coeff_ui(poly0, 1, 1);

    bool r0 = nmod_poly_is_irreducible(poly0);
    bool r1 = is_irreducible_berlekamp(poly1);
    bool r2 = is_irreducible_rabin(poly1);
    bool r3 = is_irreducible_benor(poly1);

    std::vector<ulong> gear = {1, 1};
    while (gear.size() <= 8) {
        if (r1 != r0 || r2 != r0 || r3 != r0) {
            std::cout << poly1
                      << ", expected: " << r0
                      << ", berlekamp: " << r1
                      << ", rabin: " << r2
                      << ", benor: " << r3
                      << std::endl;
        }

        gear[0] = (gear[0] + 1) % field->base();
        bool wrapped = !gear[0];
        for (uintmax_t j = 1; wrapped && j < gear.size(); ++j) {
            gear[j] = (gear[j] + 1) % field->base();
            wrapped = !gear[j];
        }
        if (wrapped) {
            gear.push_back(1);
        }

        for (slong j = 0; j < gear.size(); ++j) {
            nmod_poly_set_coeff_ui(poly0, j, gear[j]);
        }
        poly1 = gear;

        r0 = nmod_poly_is_irreducible(poly0);
        r1 = is_irreducible_berlekamp(poly1);
        r2 = is_irreducible_rabin(poly1);
        r3 = is_irreducible_benor(poly1);
    }
    return 0;
}
