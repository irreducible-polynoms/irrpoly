#include <irrpoly.h>

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

using namespace irrpoly;

[[nodiscard]]
auto fill_data(gf field, uintmax_t n) -> std::vector<gfpoly> {
    n += 2 * std::thread::hardware_concurrency();

    std::vector<gfpoly> data;
    data.push_back(gfpoly(field, {1}));
    data.push_back(gfpoly(field, {0, 1}));

    auto input = [&]() -> gfpoly {
        static std::vector<uintmax_t> gear = {0};
        static bool wrapped = false;
        gear[0] = (gear[0] + 1) % field->base();
        wrapped = !gear[0];
        for (uintmax_t i = 1; wrapped && i < gear.size(); ++i) {
            gear[i] = (gear[i] + 1) % field->base();
            wrapped = !gear[i];
        }
        if (wrapped) {
            gear.push_back(0);
        }
        std::vector<uintmax_t> copy(gear);
        copy.push_back(1);
        return gfpoly(field, copy);
    };

    while (n > 0) {
        data.emplace_back(input());
        if (is_primitive(data.back())) {
            --n;
        }
    }

    return data;
}

void bench_berlekamp(const std::vector<gfpoly> &data, uintmax_t n) {
    for (uintmax_t i = 0; n > 0; ++i) {
        if (is_irreducible_berlekamp(data[i])) {
            --n;
        }
    }
}

void bench_rabin(const std::vector<gfpoly> &data, uintmax_t n) {
    for (uintmax_t i = 0; n > 0; ++i) {
        if (is_irreducible_rabin(data[i])) {
            --n;
        }
    }
}

void bench_benor(const std::vector<gfpoly> &data, uintmax_t n) {
    for (uintmax_t i = 0; n > 0; ++i) {
        if (is_irreducible_benor(data[i])) {
            --n;
        }
    }
}

void bench_definition(const std::vector<gfpoly> &data, uintmax_t n) {
    for (uintmax_t i = 0; n > 0; ++i) {
        if (is_primitive_definition(data[i])) {
            --n;
        }
    }
}

void bench_irreducible(const std::vector<gfpoly> &data, uintmax_t n) {
    for (uintmax_t i = 0; n > 0; ++i) {
        if (is_irreducible(data[i])) {
            --n;
        }
    }
}

void bench_primitive(const std::vector<gfpoly> &data, uintmax_t n) {
    for (uintmax_t i = 0; n > 0; ++i) {
        if (is_primitive(data[i])) {
            --n;
        }
    }
}

static multithread::polychecker ch;

void bench_irr_multithread(const std::vector<gfpoly> &data, uintmax_t n) {
    uintmax_t i = 0;
    auto input = [&]() -> gfpoly {
        return data[i++];
    };

    auto check = multithread::make_check_func(
        multithread::irreducible_method::recommended,
        multithread::primitive_method::nil);

    auto callback = [&](const gfpoly &/*poly*/,
                        const typename multithread::check_result &result)
        -> bool {
        if (result.irreducible) {
            --n;
        }
        return !n;
    };

    ch.chain(input, check, callback);
}

void bench_prim_multithread(const std::vector<gfpoly> &data, uintmax_t n) {
    uintmax_t i = 0;
    auto input = [&]() -> gfpoly {
        return data[i++];
    };

    auto check = multithread::make_check_func(
        multithread::irreducible_method::nil,
        multithread::primitive_method::recommended);

    auto callback = [&](const gfpoly &/*poly*/,
                        const typename multithread::check_result &result)
        -> bool {
        if (result.primitive) {
            --n;
        }
        return !n;
    };

    ch.chain(input, check, callback);
}

/**
 * BENCHMARK functionality of Catch2 is currently under development, so there are some problems.
 * To see the correct result comment out all SECTIONs except one and get the result for that one.
 * Then repeat for others.
 */
TEST_CASE("speed test") {
    SECTION("gf2 200") {
        uintmax_t P = 2, N = 200;
        auto data = fill_data(std::move(make_gf(P)), N);
        BENCHMARK("gf2 200 berlekamp") { bench_berlekamp(data, N); };
        BENCHMARK("gf2 200 rabin") { bench_rabin(data, N); };
        BENCHMARK("gf2 200 benor")  { bench_benor(data, N); };
        BENCHMARK("gf2 200 primitive") { bench_definition(data, N); };
        BENCHMARK("gf2 200 recommended_irreducible") {
            bench_irreducible(data, N);
        };
        BENCHMARK("gf2 200 recommended_primitive") {
            bench_primitive(data, N);
        };
        BENCHMARK("gf2 200 multithread_irreducible") {
            bench_irr_multithread(data, N);
        };
        BENCHMARK("gf2 200 multithread_primitive") {
            bench_prim_multithread(data, N);
        };
    }
    SECTION("gf3 300") {
        uintmax_t P = 3, N = 300;
        auto data = fill_data(std::move(make_gf(P)), N);
        BENCHMARK("gf3 300 berlekamp") { bench_berlekamp(data, N); };
        BENCHMARK("gf3 300 rabin") { bench_rabin(data, N); };
        BENCHMARK("gf3 300 benor")  { bench_benor(data, N); };
        BENCHMARK("gf3 300 primitive") { bench_definition(data, N); };
        BENCHMARK("gf3 300 recommended_irreducible") {
            bench_irreducible(data, N);
        };
        BENCHMARK("gf3 300 recommended_primitive") {
            bench_primitive(data, N);
        };
        BENCHMARK("gf3 300 multithread_irreducible") {
            bench_irr_multithread(data, N);
        };
        BENCHMARK("gf3 300 multithread_primitive") {
            bench_prim_multithread(data, N);
        };
    }
    SECTION("gf5 400") {
        uintmax_t P = 5, N = 400;
        auto data = fill_data(std::move(make_gf(P)), N);
        BENCHMARK("gf5 400 berlekamp") { bench_berlekamp(data, N); };
        BENCHMARK("gf5 400 rabin") { bench_rabin(data, N); };
        BENCHMARK("gf5 400 benor")  { bench_benor(data, N); };
        BENCHMARK("gf5 400 primitive") { bench_definition(data, N); };
        BENCHMARK("gf5 400 recommended_irreducible") {
            bench_irreducible(data, N);
        };
        BENCHMARK("gf5 400 recommended_primitive") {
            bench_primitive(data, N);
        };
        BENCHMARK("gf5 400 multithread_irreducible") {
            bench_irr_multithread(data, N);
        };
        BENCHMARK("gf5 400 multithread_primitive") {
            bench_prim_multithread(data, N);
        };
    }
    SECTION("gf7 500") {
        uintmax_t P = 7, N = 500;
        auto data = fill_data(std::move(make_gf(P)), N);
        BENCHMARK("gf7 500 berlekamp") { bench_berlekamp(data, N); };
        BENCHMARK("gf7 500 rabin") { bench_rabin(data, N); };
        BENCHMARK("gf7 500 benor")  { bench_benor(data, N); };
        BENCHMARK("gf7 500 primitive") { bench_definition(data, N); };
        BENCHMARK("gf7 500 recommended_irreducible") {
            bench_irreducible(data, N);
        };
        BENCHMARK("gf7 500 recommended_primitive") {
            bench_primitive(data, N);
        };
        BENCHMARK("gf7 500 multithread_irreducible") {
            bench_irr_multithread(data, N);
        };
        BENCHMARK("gf7 500 multithread_primitive") {
            bench_prim_multithread(data, N);
        };
    }
    SECTION("gf11 600") {
        uintmax_t P = 11, N = 600;
        auto data = fill_data(std::move(make_gf(P)), N);
        BENCHMARK("gf11 600 berlekamp") { bench_berlekamp(data, N); };
        BENCHMARK("gf11 600 rabin") { bench_rabin(data, N); };
        BENCHMARK("gf11 600 benor")  { bench_benor(data, N); };
        BENCHMARK("gf11 600 primitive") { bench_definition(data, N); };
        BENCHMARK("gf11 600 recommended_irreducible") {
            bench_irreducible(data, N);
        };
        BENCHMARK("gf11 600 recommended_primitive") {
            bench_primitive(data, N);
        };
        BENCHMARK("gf11 600 multithread_irreducible") {
            bench_irr_multithread(data, N);
        };
        BENCHMARK("gf11 600 multithread_primitive") {
            bench_prim_multithread(data, N);
        };
    }
}
