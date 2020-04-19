#include <irrpoly/gfcheck.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

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
