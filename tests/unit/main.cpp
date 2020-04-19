#include <irrpoly.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <utility>

using namespace irrpoly;

TEST_CASE("gf could be constructed when it exists", "[gf]") {
    SECTION("empty field") {
        REQUIRE_THROWS(make_gf(0));
    }SECTION("field with only zero") {
        REQUIRE_THROWS(make_gf(1));
    }SECTION("existing field") {
        REQUIRE_NOTHROW(make_gf(2));
    }SECTION("non existing field") {
        REQUIRE_THROWS(make_gf(4));
    }SECTION("too large field") {
        REQUIRE_THROWS(make_gf(INTMAX_MAX));
    }
}

TEST_CASE("gf comparison works", "[gf]") {
    auto gf2 = make_gf(2);
    SECTION("equal to self") {
        REQUIRE(gf2 == gf2);
        REQUIRE_FALSE(gf2 != gf2);
    }SECTION("equal only to same") {
        auto same = make_gf(2);
        REQUIRE(gf2 == same);
        REQUIRE_FALSE(gf2 != same);
    }SECTION("not equal with other") {
        auto other = make_gf(3);
        REQUIRE_FALSE(gf2 == other);
        REQUIRE(gf2 != other);
    }
}

TEST_CASE("gf methods work", "[gf]") {
    auto gf5 = make_gf(5);
    SECTION("method base() works") {
        REQUIRE(gf5->base() == 5);
    }SECTION("method mul_inv() works") {
        REQUIRE_THROWS(gf5->mul_inv(0));
        REQUIRE(gf5->mul_inv(1) == 1);
        REQUIRE(gf5->mul_inv(2) == 3);
        REQUIRE(gf5->mul_inv(3) == 2);
        REQUIRE(gf5->mul_inv(4) == 4);
    }
}

TEST_CASE("gfn could be created correctly", "[gfn]") {
    auto gf5 = make_gf(5);
    SECTION("field remains the same") {
        SECTION("for direct constructed") {
            REQUIRE(gfn(gf5).field() == gf5);
            REQUIRE(gfn(gf5, 3).field() == gf5);
        }SECTION("for copy constructed") {
            auto num = gfn(gf5);
            REQUIRE(gfn(num).field() == gf5);
            REQUIRE(gfn(gfn(gf5)).field() == gf5);
        }SECTION("for randomly picked") {
            REQUIRE(gfn::random(gf5).field() == gf5);
        }SECTION("for single value") {
            REQUIRE(make_gfn(gf5, 3).field() == gf5);
        }SECTION("for initializer list") {
            auto vec = make_gfn(gf5, {3, 4, 5});
            REQUIRE(vec.size() == 3);
            REQUIRE(vec[0].field() == gf5);
            REQUIRE(vec[1].field() == gf5);
            REQUIRE(vec[2].field() == gf5);
        }SECTION("for vector") {
            auto vec = make_gfn(gf5, std::vector<uintmax_t>({3, 4, 5}));
            REQUIRE(vec.size() == 3);
            REQUIRE(vec[0].field() == gf5);
            REQUIRE(vec[1].field() == gf5);
            REQUIRE(vec[2].field() == gf5);
        }
    }SECTION("value is normalized") {
        SECTION("for direct constructed") {
            REQUIRE(gfn(gf5).data() == 0);
            REQUIRE(gfn(gf5, 7).data() == 2);
        }SECTION("for copy constructed") {
            auto num = gfn(gf5, 8);
            REQUIRE(gfn(num).data() == 3);
            REQUIRE(gfn(gfn(gf5, 9)).data() == 4);
        }SECTION("for randomly picked") {
            for (auto i : {0,1,2,3,4}) {
                REQUIRE(gfn::random(gf5).data() < 5);
            }
        }SECTION("for single value") {
            REQUIRE(make_gfn(gf5, 6).data() == 1);
        }SECTION("for initializer list") {
            auto vec = make_gfn(gf5, {6, 7, 8});
            REQUIRE(vec.size() == 3);
            REQUIRE(vec[0].data() == 1);
            REQUIRE(vec[1].data() == 2);
            REQUIRE(vec[2].data() == 3);
        }SECTION("for vector") {
            auto vec = make_gfn(gf5, std::vector<uintmax_t>({6, 7, 8}));
            REQUIRE(vec.size() == 3);
            REQUIRE(vec[0].data() == 1);
            REQUIRE(vec[1].data() == 2);
            REQUIRE(vec[2].data() == 3);
        }
    }
}
