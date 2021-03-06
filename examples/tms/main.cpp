// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include <irrpoly.h>

#include <iostream>
#include <algorithm>

using namespace irrpoly;

/// This function generates the sequence of irreducible polynomials over GF[2]
/// of growing degree, required sequence length is passed as argument.
[[nodiscard]]
auto generate_irreducible(uintmax_t num) -> std::vector<gfpoly> {
    auto gf2 = make_gf(2);

    std::vector<gfpoly> res;
    if (num == 0) {
        return res;
    }
    res.reserve(num);
    res.emplace_back(gfpoly(gf2, {0, 1}));
    if (num == 1) {
        return res;
    }

    multithread::polychecker ch;

    auto input = [&]() -> gfpoly {
        static uintmax_t index = 1;
        uint8_t degree = 0;
        for (uint8_t i = 1; index >> i; ++i, ++degree);
        std::vector<uintmax_t> res;
        res.reserve(degree + 2);
        res.emplace_back(1);
        for (uint8_t i = 0; i <= degree; ++i) {
            res.emplace_back((index & (1ull << i)) ? 1 : 0);
        }
        ++index;
        return gfpoly(gf2, res);
    };

    auto check = multithread::make_check_func(
        multithread::irreducible_method::berlekamp,
        multithread::primitive_method::nil);

    uintmax_t n = num - 1;
    auto callback = [&](const gfpoly &poly,
                        const typename multithread::check_result &result)
        -> bool {
        if (result.irreducible) {
            --n;
            res.emplace_back(poly);
        }
        return !n;
    };

    // The last argument `false` says that we want to collect all results of checks started
    // by default this argument is equal to `true` and after `callback` returns `true`
    // all results are discarded. Here we need all results for proper sequence.
    ch.chain(input, check, callback, false);

    std::sort(res.begin(), res.end(),
              [](const gfpoly &a, const gfpoly &b) {
                  if (a.degree() == b.degree()) {
                      for (auto i = a.degree(); i > 0; --i) {
                          if (a[i] != b[i]) {
                              return a[i] < b[i];
                          }
                      }
                      return a[0] < b[0];
                  }
                  return a.degree() < b.degree();
              });

    while (res.size() > num) {
        res.pop_back();
    }

    return res;
}

auto main() -> int {
    auto poly = generate_irreducible(5);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }
    return 0;
}
