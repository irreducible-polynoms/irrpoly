#include <irrpoly.h>

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <algorithm>

using namespace irrpoly;

static auto gf2 = make_gf(2);
static multithread::polychecker ch;

[[nodiscard]]
std::vector<gfpoly> generate_irreducible(uint64_t num, multithread::irreducible_method method) {
    std::vector<gfpoly> res;
    if (num == 0) {
        return res;
    }
    res.reserve(num);
    res.emplace_back(gfpoly(gf2, {0, 1}));
    if (num == 1) {
        return res;
    }

    auto input = [&]() -> gfpoly {
        static uint64_t index = 1;
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
        method,
        multithread::primitive_method::nil);

    auto callback = [&](const gfpoly &poly, const typename multithread::result_type &result) -> bool {
        if (!result.irreducible) {
            return false;
        }
        res.emplace_back(poly);
        return true;
    };

    ch.check(num - 1, input, check, callback, false);

    std::sort(res.begin(), res.end(), [](const gfpoly &a, const gfpoly &b) {
        if (a.degree() == b.degree()) {
            for (uint64_t i = a.degree(); i > 0; --i) {
                if (a[i] != b[i]) {
                    return a[i] < b[i];
                }
            }
            return a[0] < b[0];
        } else {
            return a.degree() < b.degree();
        }
    });
    while (res.size() > num) {
        res.pop_back();
    }

    return res;
}

[[nodiscard]]
std::vector<gfpoly> generate_primitive(uint64_t num) {
    std::vector<gfpoly> res;
    if (num == 0) {
        return res;
    }
    res.reserve(num);
    res.emplace_back(gfpoly(gf2, {0, 1}));
    if (num == 1) {
        return res;
    }

    auto input = [&]() -> gfpoly {
        static uint64_t index = 1;
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
        multithread::irreducible_method::nil,
        multithread::primitive_method::definition);

    auto callback = [&](const gfpoly &poly, const typename multithread::result_type &result) -> bool {
        if (!result.primitive) {
            return false;
        }
        res.emplace_back(poly);
        return true;
    };

    ch.check(num - 1, input, check, callback, false);

    std::sort(res.begin(), res.end(), [](const gfpoly &a, const gfpoly &b) {
        if (a.degree() == b.degree()) {
            for (uint64_t i = a.degree(); i > 0; --i) {
                if (a[i] != b[i]) {
                    return a[i] < b[i];
                }
            }
            return a[0] < b[0];
        } else {
            return a.degree() < b.degree();
        }
    });
    while (res.size() > num) {
        res.pop_back();
    }

    return res;
}

TEST_CASE("speed test") {
    BENCHMARK("berlekamp", i) {
                                  return generate_irreducible(i, multithread::irreducible_method::berlekamp);
                              };
    BENCHMARK("rabin", i) {
                              return generate_irreducible(i, multithread::irreducible_method::rabin);
                          };
    BENCHMARK("benor", i) {
                              return generate_irreducible(i, multithread::irreducible_method::benor);
                          };
    BENCHMARK("primitive", i) { return generate_primitive(i); };
}
