#include <irrpoly.h>

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

TEST_CASE("gfn could be created correctly", "[gfn]") {
    auto gf5 = make_gf(5);
    SECTION("for direct constructed") {
        SECTION("field remains the same") {
            REQUIRE(gfn(gf5).field() == gf5);
            REQUIRE(gfn(gf5, 3).field() == gf5);
        }SECTION("value is normalized") {
            REQUIRE(gfn(gf5).value() == 0);
            REQUIRE(gfn(gf5, 7).value() == 2);
        }
    }SECTION("for randomly picked") {
        SECTION("field remains the same") {
            REQUIRE(gfn::random(gf5).field() == gf5);
        }SECTION("value is normalized") {
            for (auto _ : {0, 1, 2, 3, 4}) {
                REQUIRE(gfn::random(gf5).value() < 5);
            }
        }
    }SECTION("for copied") {
        auto num = gfn(gf5, 2);
        SECTION("field remains the same") {
            auto field_before = num.field();
            num = gfn(gf5, 3);
            REQUIRE(field_before == num.field());
        }SECTION("value is normalized") {
            num = 10;
            REQUIRE(num.value() < 5);
        }
    }
}

TEST_CASE("gfn comparison works", "[gfn]") {
    auto gf5 = make_gf(5);
    auto num = gfn(gf5, 2);
    SECTION("comparison with numbers") {
        REQUIRE(num < gfn(gf5, 8));
        REQUIRE(num < 8);
        REQUIRE(6 < num);
    }SECTION("comparison with zero") {
        REQUIRE(!num.is_zero());
        REQUIRE(num.set_zero().is_zero());
    }SECTION("bool coerce") {
        REQUIRE(num);
        REQUIRE(!num.set_zero());
    }
}

TEST_CASE("gfn operations work", "[gfn]") {
    auto gf5 = make_gf(5);
    SECTION("sum works") {
        REQUIRE(gfn(gf5, 2) + gfn(gf5, 3) == 0);
        REQUIRE(2 + gfn(gf5, 3) == 0);
        REQUIRE(gfn(gf5, 2) + 3 == 0);
        REQUIRE(++gfn(gf5, 2) == 3);
        auto num = gfn(gf5, 2);
        REQUIRE(num++ == 2);
        REQUIRE(num == 3);
        num += 4;
        REQUIRE(num == 2);
        num += gfn(gf5, 2);
        REQUIRE(num == 4);
        REQUIRE(+num == 4);
    }SECTION("sub works") {
        REQUIRE(gfn(gf5, 2) - gfn(gf5, 3) == 4);
        REQUIRE(2 - gfn(gf5, 3) == 4);
        REQUIRE(gfn(gf5, 2) - 3 == 4);
        REQUIRE(--gfn(gf5, 2) == 1);
        auto num = gfn(gf5, 2);
        REQUIRE(num-- == 2);
        REQUIRE(num == 1);
        num -= 4;
        REQUIRE(num == 2);
        num -= gfn(gf5, 3);
        REQUIRE(num == 4);
        REQUIRE(-num == 1);
    }SECTION("mul works") {
        REQUIRE(gfn(gf5, 2) * gfn(gf5, 3) == 1);
        REQUIRE(2 * gfn(gf5, 3) == 1);
        REQUIRE(gfn(gf5, 2) * 3 == 1);
        auto num = gfn(gf5, 2);
        num *= 4;
        REQUIRE(num == 3);
        num *= gfn(gf5, 2);
        REQUIRE(num == 1);
    }SECTION("div works") {
        REQUIRE(gfn(gf5, 2) / gfn(gf5, 3) == 4);
        REQUIRE(2 / gfn(gf5, 3) == 4);
        REQUIRE(gfn(gf5, 2) / 3 == 4);
        auto num = gfn(gf5, 2);
        num /= 4;
        REQUIRE(num == 3);
        num /= gfn(gf5, 2);
        REQUIRE(num == 4);
        REQUIRE_THROWS(num / 0);
    }
}

TEST_CASE("gfpoly could be constructed correctly", "[gfpoly]") {
    auto gf5 = make_gf(5);
    std::vector<uintmax_t> etalon = {0, 1, 2, 3, 4, 0, 1};
    SECTION("empty") {
        auto poly = gfpoly(gf5);
        REQUIRE(poly.value().empty());
        REQUIRE(poly.size() == 0);
        REQUIRE_THROWS(poly.degree());
        REQUIRE(poly.field() == gf5);
        REQUIRE(poly.base() == 5);
    }SECTION("from initializer list") {
        auto poly = gfpoly(gf5, {0, 1, 2, 3, 4, 5, 6});
        REQUIRE(poly.value() == etalon);
        REQUIRE(poly.size() == etalon.size());
        REQUIRE(poly.degree() == etalon.size() - 1);
        REQUIRE(poly.field() == gf5);
        REQUIRE(poly.base() == 5);
    }SECTION("from vector") {
        std::vector<uintmax_t> vec = {0, 1, 2, 3, 4, 5, 6};
        auto poly = gfpoly(gf5, vec);
        REQUIRE(poly.value() == etalon);
        REQUIRE(poly.size() == etalon.size());
        REQUIRE(poly.degree() == etalon.size() - 1);
        REQUIRE(poly.field() == gf5);
        REQUIRE(poly.base() == 5);
    }SECTION("from number") {
        auto poly = gfpoly(gf5, 7);
        REQUIRE(poly.size() == 1);
        REQUIRE(poly.value()[0] == 2);
        REQUIRE(poly.degree() == 0);
        REQUIRE(poly.field() == gf5);
        REQUIRE(poly.base() == 5);
    }SECTION("from gfn") {
        auto poly = gfpoly(gfn(gf5, 7));
        REQUIRE(poly.size() == 1);
        REQUIRE(poly.value()[0] == 2);
        REQUIRE(poly.degree() == 0);
        REQUIRE(poly.field() == gf5);
        REQUIRE(poly.base() == 5);
    }SECTION("random") {
        for (auto i : {0, 1, 2, 3, 4}) {
            auto poly = gfpoly::random(gf5, i);
            auto normalized = gfpoly(poly).normalize();
            REQUIRE(poly.value() == normalized.value());
            REQUIRE(poly.size() == i + 1);
            REQUIRE(poly.degree() == i);
            REQUIRE(poly.field() == gf5);
            REQUIRE(poly.base() == 5);
        }
    }
}

TEST_CASE("gfpoly zero comparison works", "[gfpoly]") {
    auto gf5 = make_gf(5);
    auto poly = gfpoly::random(gf5, 2);
    REQUIRE(!poly.is_zero());
    REQUIRE(poly);
    poly.set_zero();
    REQUIRE(poly.is_zero());
    REQUIRE(!poly);
    REQUIRE(poly.value().empty());
}

TEST_CASE("gfpoly input works correctly", "[gfpoly]") {
    auto gf5 = make_gf(5);
    auto poly = gfpoly(gf5);
    REQUIRE(std::stringstream("{0, 1, 2 3, 4, 5, 6\n} ") >> poly);
    REQUIRE(poly.value() == std::vector<uintmax_t>({0, 1, 2, 3, 4, 0, 1}));
    REQUIRE_THROWS(std::stringstream("{0, 1, ") >> poly);
    REQUIRE_THROWS(std::stringstream("0, 1}") >> poly);
    REQUIRE_THROWS(std::stringstream("{-0, 1}") >> poly);
}

TEST_CASE("gfpoly operations work correctly", "[gfpoly]") {
    auto gf5 = make_gf(5);
    auto poly = gfpoly(gf5, {0, 1, 2, 3, 4});
    SECTION("rs works") {
        auto p = poly >> 2;
        REQUIRE(p == gfpoly(gf5, {2, 3, 4}));
        p >>= 2;
        REQUIRE(p == gfpoly(gf5, 4));
    }SECTION("ls works") {
        auto p = poly << 1;
        REQUIRE(p == gfpoly(gf5, {0, 0, 1, 2, 3, 4}));
        p <<= 1;
        REQUIRE(p == gfpoly(gf5, {0, 0, 0, 1, 2, 3, 4}));
    }SECTION("add works") {
        REQUIRE(poly + gfpoly(gf5, {1, 2, 3, 3, 2, 1})
                    == gfpoly(gf5, {1, 3, 0, 1, 1, 1}));
        poly += gfpoly(gf5, {1, 2, 3, 3, 2, 1});
        REQUIRE(poly == gfpoly(gf5, {1, 3, 0, 1, 1, 1}));
        REQUIRE(poly + 2 == gfpoly(gf5, {3, 3, 0, 1, 1, 1}));
        REQUIRE(2 + poly == gfpoly(gf5, {3, 3, 0, 1, 1, 1}));
        poly += 2;
        REQUIRE(poly == gfpoly(gf5, {3, 3, 0, 1, 1, 1}));
        REQUIRE(poly + gfn(gf5, 2) == gfpoly(gf5, {0, 3, 0, 1, 1, 1}));
        REQUIRE(gfn(gf5, 2) + poly == gfpoly(gf5, {0, 3, 0, 1, 1, 1}));
        poly += gfn(gf5, 2);
        REQUIRE(poly == gfpoly(gf5, {0, 3, 0, 1, 1, 1}));
    }SECTION("sub works") {
        REQUIRE(-poly == gfpoly(gf5, {0, 4, 3, 2, 1}));
        REQUIRE(poly - gfpoly(gf5, {1, 2, 3, 3, 2, 1})
                    == gfpoly(gf5, {4, 4, 4, 0, 2, 4}));
        poly -= gfpoly(gf5, {1, 2, 3, 3, 2, 1});
        REQUIRE(poly == gfpoly(gf5, {4, 4, 4, 0, 2, 4}));
        REQUIRE(poly - 2 == gfpoly(gf5, {2, 4, 4, 0, 2, 4}));
        REQUIRE(2 - poly == gfpoly(gf5, {3, 1, 1, 0, 3, 1}));
        poly -= 2;
        REQUIRE(poly == gfpoly(gf5, {2, 4, 4, 0, 2, 4}));
        REQUIRE(poly - gfn(gf5, 2) == gfpoly(gf5, {0, 4, 4, 0, 2, 4}));
        REQUIRE(gfn(gf5, 2) - poly == gfpoly(gf5, {0, 1, 1, 0, 3, 1}));
        poly -= gfn(gf5, 2);
        REQUIRE(poly == gfpoly(gf5, {0, 4, 4, 0, 2, 4}));
    }SECTION("mul works") {
        REQUIRE(poly * gfpoly(gf5, {1, 2})
                    == gfpoly(gf5, {0, 1, 4, 2, 0, 3}));
        poly *= gfpoly(gf5, {1, 2});
        REQUIRE(poly == gfpoly(gf5, {0, 1, 4, 2, 0, 3}));
        REQUIRE(poly * 2 == gfpoly(gf5, {0, 2, 3, 4, 0, 1}));
        REQUIRE(2 * poly == gfpoly(gf5, {0, 2, 3, 4, 0, 1}));
        poly *= 2;
        REQUIRE(poly == gfpoly(gf5, {0, 2, 3, 4, 0, 1}));
        REQUIRE(poly * gfn(gf5, 2) == gfpoly(gf5, {0, 4, 1, 3, 0, 2}));
        REQUIRE(gfn(gf5, 2) * poly == gfpoly(gf5, {0, 4, 1, 3, 0, 2}));
        poly *= gfn(gf5, 2);
        REQUIRE(poly == gfpoly(gf5, {0, 4, 1, 3, 0, 2}));
    }SECTION("div works") {
        REQUIRE(poly / gfpoly(gf5, {1, 1, 1})
                    == gfpoly(gf5, {4, 4, 4}));
        poly /= gfpoly(gf5, {1, 1, 1});
        REQUIRE(poly == gfpoly(gf5, {4, 4, 4}));
        REQUIRE(poly / 2 == gfpoly(gf5, {2, 2, 2}));
        poly /= 2;
        REQUIRE(poly == gfpoly(gf5, {2, 2, 2}));
        REQUIRE(poly / gfn(gf5, 2) == gfpoly(gf5, {1, 1, 1}));
        poly /= gfn(gf5, 2);
        REQUIRE(poly == gfpoly(gf5, {1, 1, 1}));
    }SECTION("rem works") {
        REQUIRE(poly % gfpoly(gf5, {1, 1, 1})
                    == gfpoly(gf5, {1, 3}));
        poly %= gfpoly(gf5, {1, 1, 1});
        REQUIRE(poly == gfpoly(gf5, {1, 3}));
        REQUIRE(poly % 2 == gfpoly(gf5));
        poly %= 2;
        REQUIRE(poly == gfpoly(gf5));
        poly = gfpoly(gf5, {1, 3});
        REQUIRE(poly % gfn(gf5, 2) == gfpoly(gf5));
        poly %= gfn(gf5, 2);
        REQUIRE(poly == gfpoly(gf5));
    }
}
