// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#define IRRPOLY_RELEASE_CHECKED
#include <irrpoly.h>

#include <iostream>
#include <thread>

using namespace irrpoly;

/// Here is the template of function for using the multithread pipeline.
[[nodiscard]]
auto generate_irreducible(
    const uintmax_t base,
    uintmax_t num,
    const uintmax_t degree,
    const typename multithread::irreducible_method irr_meth,
    const typename multithread::primitive_method prim_meth,
    const unsigned threads_num
) -> std::vector<gfpoly> {
    std::vector<gfpoly> arr;
    arr.reserve(num);

    multithread::polychecker ch(threads_num);

    auto field = make_gf(base);
    auto input = [&]() -> gfpoly {
        return gfpoly::random(field, degree);
    };

    auto check = multithread::make_check_func(irr_meth, prim_meth);

    auto callback = [&](const gfpoly &poly, const typename multithread::check_result &res) -> bool {
        if (res.irreducible) {
            --num;
            arr.emplace_back(poly);
        }
        return !num;
    };

    ch.chain(input, check, callback);

    return arr;
}

auto main() -> int {
    const uintmax_t base = 2; //< Galois field base
    const uintmax_t num = 3; //< number of polynomials to find
    const uintmax_t degree = 5; // degree of polynomials to find
    const auto irr_meth = multithread::irreducible_method::benor; // irreducibility test to use
    const auto prim_meth = multithread::primitive_method::nil; // primitivity test to use
    const unsigned threads_num = std::thread::hardware_concurrency(); // number of threads to use

    auto poly = generate_irreducible(base, num, degree, irr_meth, prim_meth, threads_num);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }

    return 0;
}
