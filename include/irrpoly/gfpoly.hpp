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
#include <cassert>
#include <ostream>
#include <algorithm>
#include <functional>
#include <initializer_list>
#include <string>
#include <sstream>
#include <cctype>

namespace irrpoly {

/**
 * Класс gfpoly представляет многочлены над полем Галуа.
 * Основан на классе polynomial из библиотеки Boost.
 */
class gfpoly final {
private:
    gf m_field;
    mutable std::vector<gfn> m_data;
    mutable bool m_normalized;

public:
    /// Генерирует случайный многочлен над полем GF[P] заданной степени.
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

    auto value() const -> std::vector<uintmax_t> {
        normalize();
        std::vector<uintmax_t> res;
        res.reserve(size());
        for (const gfn &val : m_data) {
            res.emplace_back(val.value());
        }
        return res;
    }

    auto normalize() -> gfpoly & {
        if (!m_normalized) {
            m_data.erase(std::find_if(
                m_data.rbegin(), m_data.rend(),
                std::not_fn(std::mem_fn(&gfn::is_zero))
            ).base(), m_data.end());
            m_normalized = true;
        }
        return *this;
    }

private:
    auto normalize() const -> const gfpoly & {
        if (!m_normalized) {
            m_data.erase(std::find_if(
                m_data.rbegin(), m_data.rend(),
                std::not_fn(std::mem_fn(&gfn::is_zero))
            ).base(), m_data.end());
            m_normalized = true;
        }
        return *this;
    }

    gfpoly(const gf &field, std::vector<gfn> &&p) :
        m_field(field), m_data(std::move(p)), m_normalized(false) {
        normalize();
    }

public:
    explicit
    gfpoly(const gf &field) : m_field(field), m_data(), m_normalized(true) {}

    gfpoly(const gf &field, const std::vector<uintmax_t> &l) :
        m_field(field), m_data(), m_normalized(false) {
        m_data.reserve(l.size());
        for (uintmax_t v : l) {
            m_data.emplace_back(field, v);
        }
        normalize();
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

    gfpoly(const gfpoly &p) :
        m_field(p.m_field), m_data(p.m_data), m_normalized(p.m_normalized) {
        normalize();
    };

    gfpoly(gfpoly &&p) noexcept:
        m_field(std::move(p.m_field)), m_data(std::move(p.m_data)),
        m_normalized(p.m_normalized) {
        normalize();
    };

    auto operator=(const gfpoly &p) -> gfpoly & {
        if (this != &p) {
            assert(m_field == nullptr || m_field == p.m_field);
            m_field = p.m_field;
            m_data = p.m_data;
            m_normalized = p.m_normalized;
        }
        normalize();
        return *this;
    }

    explicit
    gfpoly(gfn value) :
        m_field(value.field()), m_data(), m_normalized(false) {
        m_data.push_back(std::move(value));
        normalize();
    }

    auto operator=(gfn value) -> gfpoly & {
        assert(m_field == nullptr || m_field == value.field());
        gfpoly copy(std::move(value));
        std::swap(*this, copy);
        return *this;
    }

    gfpoly(const gf &field, uintmax_t value) :
        m_field(field), m_data(), m_normalized(false) {
        m_data.emplace_back(field, value);
        normalize();
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
        if (!m_normalized) {
            normalize();
        }
        return m_data.size();
    }

    [[nodiscard]]
    auto degree() const -> uintmax_t {
        normalize();
        if (size() == 0) {
            throw std::logic_error("degree() is undefined for the zero polynomial");
        }
        return m_data.size() - 1;
    }

    auto operator[](uintmax_t i) -> gfn & {
        m_normalized = i + 1 == m_data.size() ? false : m_normalized;
        return m_data[i];
    }

    auto operator[](uintmax_t i) const -> const gfn & {
        return m_data[i];
    }

    [[nodiscard]]
    auto data() const -> const std::vector<gfn> & {
        normalize();
        return m_data;
    }

private:
    auto data() -> std::vector<gfn> & {
        m_normalized = false;
        return m_data;
    }

public:
    [[nodiscard]]
    auto is_zero() const -> bool {
        normalize();
        return m_data.empty();
    }

    explicit operator bool() const {
        return !is_zero();
    }

    auto set_zero() -> gfpoly & {
        m_data.clear();
        m_normalized = true;
        return *this;
    }

private:
    template<class U, class R>
    auto transform(const U &value, R op) -> gfpoly & {
        m_normalized = false;
        if (m_data.empty()) {
            m_data.resize(1, gfn(m_field));
        }
        m_data[0] = op(m_data[0], value);
        return normalize();
    }

    template<class R>
    auto transform(const gfpoly &value, R op) -> gfpoly & {
        assert(field() == value.field());
        m_normalized = false;
        if (m_data.size() < value.size()) {
            m_data.resize(value.size(), gfn(m_field));
        }
        for (uintmax_t i = 0; i < value.size(); ++i) {
            m_data[i] = op(m_data[i], value[i]);
        }
        return normalize();
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
        m_normalized = false;
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x * value; });
        return normalize();
    }

    template<class U>
    auto operator/=(const U &value) -> gfpoly & {
        m_normalized = false;
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x / value; });
        return normalize();
    }

    template<class U>
    auto operator%=(const U & /*value_type*/) -> gfpoly & {
        // We can always divide by a scalar, so there is no remainder:
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
        assert(a.field() == b.field());
        m_normalized = false;
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
        return normalize();
    }

public:
    auto operator*=(const gfpoly &value) -> gfpoly & {
        return multiply(*this, value);
    }

private:
    template<typename N>
    static void division_impl(gfpoly &q, gfpoly &u, const gfpoly &v, N n, N k) {
        assert(q.field() == u.field() && u.field() == v.field());
        q[k] = u[n + k] / v[n];
        for (N j = n + k; j > k;) {
            j--;
            u[j] -= q[k] * v[j - k];
        }
    }

    static auto division(gfpoly u, const gfpoly &v) -> std::pair<gfpoly, gfpoly> {
        assert(u.field() == v.field() && v.size() <= u.size() && v && u);
        uintmax_t const m = u.size() - 1, n = v.size() - 1;
        uintmax_t k = m - n;
        gfpoly q(u.field());
        q.data().resize(m - n + 1, gfn(q.field()));

        do {
            division_impl(q, u, v, n, k);
        } while (k-- != 0);
        u.data().resize(n, gfn(q.field())), gfn(q.field());
        return std::make_pair(q.normalize(), u.normalize());
    }

    /** Calculates a / b and a % b, returning the pair (quotient, remainder) together
      * because the same amount of computation yields both.
      * This function is not defined for division by zero: user beware.
      */
    static auto quotient_remainder(const gfpoly &dividend, const gfpoly &divisor) -> std::pair<gfpoly, gfpoly> {
        assert(dividend.field() == divisor.field() && divisor);
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

    template<typename U>
    auto operator>>=(U const &n) -> gfpoly & {
        m_data.erase(m_data.begin(), m_data.begin() + (n <= size() ? n : size()));
        return *this;
    }

    template<typename U>
    auto operator<<=(U const &n) -> gfpoly & {
        normalize();
        m_data.insert(m_data.begin(), n, gfn(m_field));
        return *this;
    }

    friend auto operator-(gfpoly /*a*/) -> gfpoly;

    friend auto operator*(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator/(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator%(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;
};

auto operator-(gfpoly a) -> gfpoly {
    a.m_normalized = false;
    std::transform(a.data().begin(), a.data().end(), a.data().begin(), std::negate());
    return a.normalize();
}

auto operator+(const gfpoly &a, const gfpoly &b) -> gfpoly {
    assert(a.field() == b.field());
    gfpoly result(a);
    result += b;
    return result;
}

auto operator+(gfpoly &&a, const gfpoly &b) -> gfpoly {
    assert(a.field() == b.field());
    a += b;
    return a;
}

auto operator+(const gfpoly &a, gfpoly &&b) -> gfpoly {
    assert(a.field() == b.field());
    b += a;
    return b;
}

auto operator+(gfpoly &&a, gfpoly &&b) -> gfpoly {
    assert(a.field() == b.field());
    a += b;
    return a;
}

auto operator-(const gfpoly &a, const gfpoly &b) -> gfpoly {
    assert(a.field() == b.field());
    gfpoly result(a);
    result -= b;
    return result;
}

auto operator-(gfpoly &&a, const gfpoly &b) -> gfpoly {
    assert(a.field() == b.field());
    a -= b;
    return a;
}

auto operator-(const gfpoly &a, gfpoly &&b) -> gfpoly {
    assert(a.field() == b.field());
    b -= a;
    return -b;
}

auto operator-(gfpoly &&a, gfpoly &&b) -> gfpoly {
    assert(a.field() == b.field());
    a -= b;
    return a;
}

auto operator*(const gfpoly &a, const gfpoly &b) -> gfpoly {
    gfpoly result(a.field());
    return result.multiply(a, b);
}

auto operator/(const gfpoly &a, const gfpoly &b) -> gfpoly {
    return gfpoly::quotient_remainder(a, b).first;
}

auto operator%(const gfpoly &a, const gfpoly &b) -> gfpoly {
    return gfpoly::quotient_remainder(a, b).second;
}

template<class U>
auto operator+(gfpoly a, const U &b) -> gfpoly {
    a += b;
    return a;
}

template<class U>
auto operator-(gfpoly a, const U &b) -> gfpoly {
    a -= b;
    return a;
}

template<class U>
auto operator*(gfpoly a, const U &b) -> gfpoly {
    a *= b;
    return a;
}

template<class U>
auto operator/(gfpoly a, const U &b) -> gfpoly {
    a /= b;
    return a;
}

template<class U>
auto operator%(const gfpoly &a, const U & /*unused*/) -> gfpoly {
    return gfpoly(a.field());
}

template<class U>
auto operator+(const U &a, gfpoly b) -> gfpoly {
    b += a;
    return b;
}

template<class U>
auto operator-(const U &a, gfpoly b) -> gfpoly {
    b -= a;
    return -b;
}

template<class U>
auto operator*(const U &a, gfpoly b) -> gfpoly {
    b *= a;
    return b;
}

auto operator==(const gfpoly &a, const gfpoly &b) -> bool {
    assert(a.field() == b.field());
    return a.data() == b.data();
}

auto operator!=(const gfpoly &a, const gfpoly &b) -> bool {
    assert(a.field() == b.field());
    return a.data() != b.data();
}

template<typename U>
auto operator>>(gfpoly a, const U &b) -> gfpoly {
    a >>= b;
    return a;
}

template<typename U>
auto operator<<(gfpoly a, const U &b) -> gfpoly {
    a <<= b;
    return a;
}

template<class charT, class traits>
auto operator<<(std::basic_ostream<charT, traits> &os, const gfpoly &poly) -> std::basic_ostream<charT, traits> & {
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
auto operator>>(std::basic_istream<charT, traits> &is, gfpoly &poly) -> std::basic_istream<charT, traits> & {
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

} // namespace irrpoly
