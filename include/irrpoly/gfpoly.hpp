/**
 * @file    gfpoly.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include "gf.hpp"

#include <vector>
#include <functional>
#include <initializer_list>
#include <string>
#include <sstream>
#include <cctype>

namespace irrpoly {

/**
 * Binary operations for two gfn instances are correctly defined only
 * when field is the same for both of them. By default this is checked
 * only in Debug configuration and no checks performed in Release to speed
 * up computations. If you are not sure in correctness of your code add
 * #define IRRPOLY_RELEASE_CHECKED before #include <irrpoly.h> to enable
 * checks for Release configuration.
 */
#if !defined(NDEBUG) || defined(IRRPOLY_RELEASE_CHECKED) // Debug or Release Checked
#define CHECK_FIELD(comparison) \
    if (!(comparison)) { \
        throw std::logic_error("field check failed"); \
    }
#else // Release
#define CHECK_FIELD(comparison)
#endif

/**
 * gfpoly represents a polynomial over Galois field.
 * This class is originally taken from Boost library but was significantly changed.
 * Get by index operation [i] returns polynomial term x^i.
 * Polynomial is either zero or reduced which means that leading term is non-zero.
 */
class gfpoly final {
private:
    gf m_field;
    std::vector<gfn> m_data; ///< polynomial coefficients

public:
    /**
     * Generates random polynomial over provided Galois field of given degree.
     */
    static
    auto random(const gf &field, uintmax_t degree) -> gfpoly {
        std::vector<gfn> data;
        data.reserve(degree + 1);
        for (uintmax_t i = 0; i < degree; ++i) {
            data.emplace_back(gfn::random(field));
        }
        data.emplace_back(field, 1);
        while (data[0].is_zero()) {
            data[0] = gfn::random(field);
        }
        return gfpoly(field, std::move(data));
    }

    /**
     * Returns copy of numeric representation of current polynomial.
     */
    [[nodiscard]]
    auto value() const -> std::vector<uintmax_t> {
        std::vector<uintmax_t> res;
        res.reserve(size());
        for (const gfn &val : m_data) {
            res.emplace_back(val.value());
        }
        return res;
    }

private:
    /**
     * Removes leading zeroes.
     */
    auto reduce() -> gfpoly & {
        m_data.erase(std::find_if(
            m_data.rbegin(), m_data.rend(),
            std::not_fn(std::mem_fn(&gfn::is_zero))
        ).base(), m_data.end());
        return *this;
    }

    /**
     * This constructor is private because it do not checks that
     * all gfn instances have the same field.
     */
    gfpoly(const gf &field, std::vector<gfn> &&p) :
        m_field(field), m_data(std::move(p)) {
        reduce();
    }

public:
    explicit
    gfpoly(const gf &field) : m_field(field), m_data() {}

    gfpoly(const gf &field, const std::vector<uintmax_t> &l) :
        m_field(field), m_data() {
        m_data.reserve(l.size());
        for (uintmax_t v : l) {
            m_data.emplace_back(field, v);
        }
        reduce();
    }

    auto operator=(const std::vector<uintmax_t> &l) -> gfpoly & {
        gfpoly copy(m_field, l);
        std::swap(*this, copy);
        return *this;
    }

    gfpoly(const gf &field, std::initializer_list<uintmax_t> l) :
        gfpoly(field, std::vector<uintmax_t>{l}) {}

    auto operator=(std::initializer_list<uintmax_t> l) -> gfpoly & {
        gfpoly copy(m_field, l);
        std::swap(*this, copy);
        return *this;
    }

    gfpoly(const gfpoly &p) = default;

    gfpoly(gfpoly &&p) = default;

    auto operator=(const gfpoly &p) -> gfpoly & {
        if (this != &p) {
            // m_field == nullptr means that gfn instance is uninitialised
            // this happens during std::move, std::swap and inside some std::vector methods
            CHECK_FIELD(m_field == nullptr || m_field == p.m_field)
            m_field = p.m_field;
            m_data = p.m_data;
        }
        return *this;
    }

    explicit
    gfpoly(gfn value) :
        m_field(value.field()), m_data() {
        if (value) {
            m_data.push_back(std::move(value));
        }
    }

    auto operator=(gfn value) -> gfpoly & {
        // m_field == nullptr means that gfn instance is uninitialised
        // this happens during std::move, std::swap and inside some std::vector methods
        CHECK_FIELD(m_field == nullptr || m_field == value.field())
        gfpoly copy(std::move(value));
        std::swap(*this, copy);
        return *this;
    }

    gfpoly(const gf &field, uintmax_t value) :
        m_field(field), m_data() {
        if (value % base() != 0) {
            m_data.emplace_back(field, value);
        }
    }

    auto operator=(uintmax_t value) -> gfpoly & {
        gfpoly copy(m_field, value);
        std::swap(*this, copy);
        return *this;
    }

    [[nodiscard]]
    auto field() const -> const gf & {
        return m_field;
    }

    [[nodiscard]]
    auto base() const -> uintmax_t {
        return m_field->base();
    }

    [[nodiscard]]
    auto size() const -> uintmax_t {
        return m_data.size();
    }

    /**
     * Returns polynomial degree. For zero polynomial degree is undefined.
     */
    [[nodiscard]]
    auto degree() const -> uintmax_t {
        if (size() == 0) {
            throw std::logic_error("degree is undefined for zero polynomial");
        }
        return m_data.size() - 1;
    }

    /**
     * Polynomial couldn't be mutated by index because it would be possible to
     * replace some term by gfn instance with field different to polynomial's one.
     */
    auto operator[](uintmax_t i) const -> const gfn & {
        return m_data[i];
    }

    [[nodiscard]]
    auto is_zero() const -> bool {
        return m_data.empty();
    }

    explicit operator bool() const {
        return !is_zero();
    }

    [[maybe_unused]]
    auto set_zero() -> gfpoly & {
        m_data.clear();
        return *this;
    }

private:
    template<class U, class R>
    auto transform(const U &value, R op) -> gfpoly & {
        if (m_data.empty()) {
            m_data.resize(1, gfn(m_field));
        }
        m_data[0] = op(m_data[0], value);
        return reduce();
    }

    template<class R>
    auto transform(const gfpoly &value, R op) -> gfpoly & {
        CHECK_FIELD(field() == value.field())
        if (m_data.size() < value.size()) {
            m_data.resize(value.size(), gfn(m_field));
        }
        for (uintmax_t i = 0; i < value.size(); ++i) {
            m_data[i] = op(m_data[i], value[i]);
        }
        return reduce();
    }

public:
    template<class U>
    auto operator+=(const U &value) -> gfpoly & {
        return transform(value, std::plus());
    }

    template<class U>
    auto operator-=(const U &value) -> gfpoly & {
        return transform(value, std::minus());
    }

    template<class U>
    auto operator*=(const U &value) -> gfpoly & {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x * value; });
        return reduce();
    }

    template<class U>
    auto operator/=(const U &value) -> gfpoly & {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x / value; });
        return reduce();
    }

    template<class U>
    auto operator%=(const U & /*value*/) -> gfpoly & {
        // we can always divide by a scalar, so there is no remainder
        return set_zero();
    }

    auto operator+=(const gfpoly &value) -> gfpoly & {
        return transform(value, std::plus());
    }

    auto operator-=(const gfpoly &value) -> gfpoly & {
        return transform(value, std::minus());
    }

private:
    auto multiply(const gfpoly &a, const gfpoly &b) -> gfpoly & {
        CHECK_FIELD(a.field() == b.field())
        if (!a || !b) {
            return set_zero();
        }
        std::vector<gfn> prod(a.size() + b.size() - 1, gfn(m_field));
        for (uintmax_t i = 0; i < a.size(); ++i) {
            for (uintmax_t j = 0; j < b.size(); ++j) {
                prod[i + j] += a.m_data[i] * b.m_data[j];
            }
        }
        m_data.swap(prod);
        return reduce();
    }

public:
    auto operator*=(const gfpoly &value) -> gfpoly & {
        return multiply(*this, value);
    }

private:
    /**
     * This is the core method of the whole class. It takes the most time during computations.
     * TODO: try implementing polynomial division using Discrete Fourier Transform.
     */
    static auto division(gfpoly u, const gfpoly &v) -> std::pair<gfpoly, gfpoly> {
        CHECK_FIELD(u.field() == v.field() && v.size() <= u.size() && v && u)
        uintmax_t const m = u.size() - 1, n = v.size() - 1;
        uintmax_t k = m - n;
        gfpoly q(u.field());
        q.m_data.resize(m - n + 1, gfn(q.field()));

        auto division_impl = [](gfpoly *q, gfpoly *u, const gfpoly &v, auto n, auto k) {
            CHECK_FIELD(q->field() == u->field() && u->field() == v.field())
            q->m_data[k] = u->m_data[n + k] / v.m_data[n];
            for (auto j = n + k; j > k;) {
                j--;
                u->m_data[j] -= q->m_data[k] * v.m_data[j - k];
            }
        };

        do {
            division_impl(&q, &u, v, n, k);
        } while (k-- != 0);
        u.m_data.resize(n, gfn(q.field())), gfn(q.field());
        return std::make_pair(q.reduce(), u.reduce());
    }

    /**
     * Calculates a / b and a % b, returning the pair (quotient, remainder) together
     * because the same amount of computation yields both.
     * This function is not defined for division by zero: user beware.
     */
    static auto quotient_remainder(const gfpoly &dividend, const gfpoly &divisor)
    -> std::pair<gfpoly, gfpoly> {
        CHECK_FIELD(dividend.field() == divisor.field() && divisor)
        if (dividend.size() < divisor.size()) {
            return std::make_pair(gfpoly(dividend.field()), dividend);
        }
        return division(dividend, divisor);
    }

public:
    auto operator/=(const gfpoly &value) -> gfpoly & {
        return *this = quotient_remainder(*this, value).first;
    }

    auto operator%=(const gfpoly &value) -> gfpoly & {
        return *this = quotient_remainder(*this, value).second;
    }

    /**
     * Logically equal to operation this /= x^n. Defined only when such devision is possible.
     */
    template<typename U>
    auto operator>>=(U const &n) -> gfpoly & {
        if (n > degree() || !std::all_of(
            m_data.begin(), m_data.begin() + n, std::mem_fn(&gfn::is_zero))) {
            throw std::logic_error("division is impossible");
        }
        m_data.erase(m_data.begin(), m_data.begin() + n);
        return *this;
    }

    /**
     * Logically equal to operation this *= x^n.
     */
    template<typename U>
    auto operator<<=(U const &n) -> gfpoly & {
        reduce();
        m_data.insert(m_data.begin(), n, gfn(m_field));
        return *this;
    }

    friend auto operator-(gfpoly /*a*/) -> gfpoly;

    friend auto operator*(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator/(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator%(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator==(const gfpoly & /*a*/, const gfpoly & /*b*/) -> bool;

    friend auto operator!=(const gfpoly & /*a*/, const gfpoly & /*b*/) -> bool;
};

inline
auto operator-(gfpoly a) -> gfpoly {
    std::transform(a.m_data.begin(), a.m_data.end(), a.m_data.begin(), std::negate());
    return a.reduce();
}

inline
auto operator+(const gfpoly &a, const gfpoly &b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    gfpoly result(a);
    result += b;
    return result;
}

inline
auto operator+(gfpoly &&a, const gfpoly &b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    a += b;
    return a;
}

inline
auto operator+(const gfpoly &a, gfpoly &&b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    b += a;
    return b;
}

inline
auto operator+(gfpoly &&a, gfpoly &&b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    a += b;
    return a;
}

inline
auto operator-(const gfpoly &a, const gfpoly &b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    gfpoly result(a);
    result -= b;
    return result;
}

inline
auto operator-(gfpoly &&a, const gfpoly &b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    a -= b;
    return a;
}

inline
auto operator-(const gfpoly &a, gfpoly &&b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    b -= a;
    return -b;
}

inline
auto operator-(gfpoly &&a, gfpoly &&b) -> gfpoly {
    CHECK_FIELD(a.field() == b.field())
    a -= b;
    return a;
}

inline
auto operator*(const gfpoly &a, const gfpoly &b) -> gfpoly {
    gfpoly result(a.field());
    return result.multiply(a, b);
}

inline
auto operator/(const gfpoly &a, const gfpoly &b) -> gfpoly {
    return gfpoly::quotient_remainder(a, b).first;
}

inline
auto operator%(const gfpoly &a, const gfpoly &b) -> gfpoly {
    return gfpoly::quotient_remainder(a, b).second;
}

template<class U>
inline
auto operator+(gfpoly a, const U &b) -> gfpoly {
    a += b;
    return a;
}

template<class U>
inline
auto operator-(gfpoly a, const U &b) -> gfpoly {
    a -= b;
    return a;
}

template<class U>
inline
auto operator*(gfpoly a, const U &b) -> gfpoly {
    a *= b;
    return a;
}

template<class U>
inline
auto operator/(gfpoly a, const U &b) -> gfpoly {
    a /= b;
    return a;
}

template<class U>
inline
auto operator%(const gfpoly &a, const U & /*unused*/) -> gfpoly {
    return gfpoly(a.field());
}

template<class U>
inline
auto operator+(const U &a, gfpoly b) -> gfpoly {
    b += a;
    return b;
}

template<class U>
inline
auto operator-(const U &a, gfpoly b) -> gfpoly {
    b -= a;
    return -b;
}

template<class U>
inline
auto operator*(const U &a, gfpoly b) -> gfpoly {
    b *= a;
    return b;
}

inline
auto operator==(const gfpoly &a, const gfpoly &b) -> bool {
    CHECK_FIELD(a.field() == b.field())
    return a.m_data == b.m_data;
}

inline
auto operator!=(const gfpoly &a, const gfpoly &b) -> bool {
    CHECK_FIELD(a.field() == b.field())
    return a.m_data != b.m_data;
}

template<typename U>
inline
auto operator>>(gfpoly a, const U &b) -> gfpoly {
    a >>= b;
    return a;
}

template<typename U>
inline
auto operator<<(gfpoly a, const U &b) -> gfpoly {
    a <<= b;
    return a;
}

template<class charT, class traits>
auto operator<<(std::basic_ostream<charT, traits> &os, const gfpoly &poly)
-> std::basic_ostream<charT, traits> & {
    os << "{ ";
    for (uintmax_t i = 0; i < poly.size(); ++i) {
        if (i) {
            os << ", ";
        }
        os << poly[i];
    }
    os << " }";
    return os;
}

template<class charT, class traits>
auto operator>>(std::basic_istream<charT, traits> &is, gfpoly &poly)
-> std::basic_istream<charT, traits> & {
    charT tmp;
    uintmax_t num = 0;
    std::string str;
    std::vector<uintmax_t> vec;
    while (is.good() && is.get(tmp) && tmp != '{') {
        if (tmp != ' ' && tmp != '\n') {
            throw std::invalid_argument("wrong input");
        }
    }
    if (tmp != '{') {
        throw std::invalid_argument("wrong input");
    }
    while (is.good() && is.get(tmp) && tmp != '}') {
        if (tmp == ',' || tmp == ' ' || tmp == '\n') {
            if (!str.empty()) {
                std::stringstream(str) >> num;
                str.clear();
                vec.emplace_back(num);
            }
        } else if (std::isdigit(tmp)) {
            str += tmp;
        } else {
            throw std::invalid_argument("wrong input");
        }
    }
    if (tmp != '}') {
        throw std::invalid_argument("wrong input");
    }
    poly = gfpoly(poly.field(), vec);
    return is;
}

#undef CHECK_FIELD

} // namespace irrpoly
