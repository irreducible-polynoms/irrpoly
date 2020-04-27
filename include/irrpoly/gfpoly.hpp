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
    gfpoly random(const gf &field, uintmax_t degree) {
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

    std::vector<uintmax_t> value() const {
        normalize();
        std::vector<uintmax_t> res;
        res.reserve(size());
        for (const gfn &val : m_data) {
            res.emplace_back(val.value());
        }
        return res;
    }

    gfpoly &normalize() {
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
    const gfpoly &normalize() const {
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

    gfpoly &operator=(const std::vector<uintmax_t> &l) {
        gfpoly copy(m_field, l);
        std::swap(*this, copy);
        return *this;
    }

    gfpoly(const gf &field, std::initializer_list<uintmax_t> l) :
        gfpoly(field, std::vector<uintmax_t>{l}) {}

    gfpoly &operator=(std::initializer_list<uintmax_t> l) {
        gfpoly copy(m_field, l);
        std::swap(*this, copy);
        return *this;
    }

    gfpoly(const gfpoly &p) :
        m_field(p.m_field), m_data(p.m_data), m_normalized(p.m_normalized) {
        normalize();
    };

    gfpoly(gfpoly &&p) noexcept :
        m_field(std::move(p.m_field)), m_data(std::move(p.m_data)),
        m_normalized(p.m_normalized) {
        normalize();
    };

    gfpoly &operator=(const gfpoly &p) {
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

    gfpoly &operator=(gfn value) {
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

    gfpoly &operator=(uintmax_t value) {
        gfpoly copy(m_field, value);
        std::swap(*this, copy);
        return *this;
    }

    [[nodiscard]]
    const gf &field() const {
        return m_field;
    }

    [[nodiscard]]
    uintmax_t base() const {
        return m_field->base();
    }

    [[nodiscard]]
    uintmax_t size() const {
        if (!m_normalized) {
            normalize();
        }
        return m_data.size();
    }

    [[nodiscard]]
    uintmax_t degree() const {
        normalize();
        if (size() == 0)
            throw std::logic_error("degree() is undefined for the zero polynomial");
        return m_data.size() - 1;
    }

    gfn &operator[](uintmax_t i) {
        m_normalized = false;
        return m_data[i];
    }

    const gfn &operator[](uintmax_t i) const {
        return m_data[i];
    }

    [[nodiscard]]
    const std::vector<gfn> &data() const {
        normalize();
        return m_data;
    }

private:
    std::vector<gfn> &data() {
        m_normalized = false;
        return m_data;
    }

public:
    [[nodiscard]]
    bool is_zero() const {
        normalize();
        return m_data.empty();
    }

    explicit operator bool() const {
        return !is_zero();
    }

    gfpoly &set_zero() {
        m_data.clear();
        m_normalized = true;
        return *this;
    }

private:
    template<class U, class R>
    gfpoly &transform(const U &value, R op) {
        m_normalized = false;
        if (m_data.empty())
            m_data.resize(1, gfn(m_field));
        m_data[0] = op(m_data[0], value);
        return normalize();
    }

    template<class R>
    gfpoly &transform(const gfpoly &value, R op) {
        assert(field() == value.field());
        m_normalized = false;
        if (m_data.size() < value.size())
            m_data.resize(value.size(), gfn(m_field));
        for (uintmax_t i = 0; i < value.size(); ++i)
            m_data[i] = op(m_data[i], value[i]);
        return normalize();
    }

public:
    template<class U>
    gfpoly &operator+=(const U &value) {
        return transform(value, std::plus());
    }

    template<class U>
    gfpoly &operator-=(const U &value) {
        return transform(value, std::minus());
    }

    template<class U>
    gfpoly &operator*=(const U &value) {
        m_normalized = false;
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x * value; });
        return normalize();
    }

    template<class U>
    gfpoly &operator/=(const U &value) {
        m_normalized = false;
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x / value; });
        return normalize();
    }

    template<class U>
    gfpoly &operator%=(const U & /*value_type*/) {
        // We can always divide by a scalar, so there is no remainder:
        return set_zero();
    }

    gfpoly &operator+=(const gfpoly &value) {
        return transform(value, std::plus());
    }

    gfpoly &operator-=(const gfpoly &value) {
        return transform(value, std::minus());
    }

private:
    gfpoly &multiply(const gfpoly &a, const gfpoly &b) {
        assert(a.field() == b.field());
        m_normalized = false;
        if (!a || !b) {
            return set_zero();
        }
        std::vector<gfn> prod(a.size() + b.size() - 1, gfn(m_field));
        for (typename std::vector<gfn>::size_type i = 0; i < a.size(); ++i)
            for (typename std::vector<gfn>::size_type j = 0; j < b.size(); ++j)
                prod[i + j] += a.m_data[i] * b.m_data[j];
        m_data.swap(prod);
        return normalize();
    }

public:
    gfpoly &operator*=(const gfpoly &value) {
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

    static std::pair<gfpoly, gfpoly> division(gfpoly u, const gfpoly &v) {
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
    static std::pair<gfpoly, gfpoly> quotient_remainder(const gfpoly &dividend, const gfpoly &divisor) {
        assert(dividend.field() == divisor.field() && divisor);
        if (dividend.size() < divisor.size())
            return std::make_pair(gfpoly(dividend.field()), dividend);
        return division(dividend, divisor);
    }

public:
    gfpoly &operator/=(const gfpoly &value) {
        return *this = quotient_remainder(*this, value).first;
    }

    gfpoly &operator%=(const gfpoly &value) {
        return *this = quotient_remainder(*this, value).second;
    }

    template<typename U>
    gfpoly &operator>>=(U const &n) {
        m_data.erase(m_data.begin(), m_data.begin() + (n <= size() ? n : size()));
        return *this;
    }

    template<typename U>
    gfpoly &operator<<=(U const &n) {
        normalize();
        m_data.insert(m_data.begin(), n, gfn(m_field));
        return *this;
    }

    friend gfpoly operator-(gfpoly);

    friend gfpoly operator*(const gfpoly &, const gfpoly &);

    friend gfpoly operator/(const gfpoly &, const gfpoly &);

    friend gfpoly operator%(const gfpoly &, const gfpoly &);
};

gfpoly operator-(gfpoly a) {
    a.m_normalized = false;
    std::transform(a.data().begin(), a.data().end(), a.data().begin(), std::negate());
    return a.normalize();
}

gfpoly operator+(const gfpoly &a, const gfpoly &b) {
    assert(a.field() == b.field());
    gfpoly result(a);
    result += b;
    return result;
}

gfpoly operator+(gfpoly &&a, const gfpoly &b) {
    assert(a.field() == b.field());
    a += b;
    return a;
}

gfpoly operator+(const gfpoly &a, gfpoly &&b) {
    assert(a.field() == b.field());
    b += a;
    return b;
}

gfpoly operator+(gfpoly &&a, gfpoly &&b) {
    assert(a.field() == b.field());
    a += b;
    return a;
}

gfpoly operator-(const gfpoly &a, const gfpoly &b) {
    assert(a.field() == b.field());
    gfpoly result(a);
    result -= b;
    return result;
}

gfpoly operator-(gfpoly &&a, const gfpoly &b) {
    assert(a.field() == b.field());
    a -= b;
    return a;
}

gfpoly operator-(const gfpoly &a, gfpoly &&b) {
    assert(a.field() == b.field());
    b -= a;
    return -b;
}

gfpoly operator-(gfpoly &&a, gfpoly &&b) {
    assert(a.field() == b.field());
    a -= b;
    return a;
}

gfpoly operator*(const gfpoly &a, const gfpoly &b) {
    gfpoly result(a.field());
    return result.multiply(a, b);
}

gfpoly operator/(const gfpoly &a, const gfpoly &b) {
    return gfpoly::quotient_remainder(a, b).first;
}

gfpoly operator%(const gfpoly &a, const gfpoly &b) {
    return gfpoly::quotient_remainder(a, b).second;
}

template<class U>
gfpoly operator+(gfpoly a, const U &b) {
    a += b;
    return a;
}

template<class U>
gfpoly operator-(gfpoly a, const U &b) {
    a -= b;
    return a;
}

template<class U>
gfpoly operator*(gfpoly a, const U &b) {
    a *= b;
    return a;
}

template<class U>
gfpoly operator/(gfpoly a, const U &b) {
    a /= b;
    return a;
}

template<class U>
gfpoly operator%(const gfpoly &a, const U &) {
    return gfpoly(a.field());
}

template<class U>
gfpoly operator+(const U &a, gfpoly b) {
    b += a;
    return b;
}

template<class U>
gfpoly operator-(const U &a, gfpoly b) {
    b -= a;
    return -b;
}

template<class U>
gfpoly operator*(const U &a, gfpoly b) {
    b *= a;
    return b;
}

bool operator==(const gfpoly &a, const gfpoly &b) {
    assert(a.field() == b.field());
    return a.data() == b.data();
}

bool operator!=(const gfpoly &a, const gfpoly &b) {
    assert(a.field() == b.field());
    return a.data() != b.data();
}

template<typename U>
gfpoly operator>>(gfpoly a, const U &b) {
    a >>= b;
    return a;
}

template<typename U>
gfpoly operator<<(gfpoly a, const U &b) {
    a <<= b;
    return a;
}

template<class charT, class traits>
std::basic_ostream<charT, traits> &
operator<<(std::basic_ostream<charT, traits> &os, const gfpoly &poly) {
    os << "{ ";
    for (uintmax_t i = 0; i < poly.size(); ++i) {
        if (i)
            os << ", ";
        os << poly[i];
    }
    os << " }";
    return os;
}

template<class charT, class traits>
std::basic_istream<charT, traits> &
operator>>(std::basic_istream<charT, traits> &is, gfpoly &poly) {
    charT tmp;
    uintmax_t num;
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
