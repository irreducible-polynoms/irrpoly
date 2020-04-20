/**
 * @file    gfpoly.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#ifndef GFPOLY_HPP
#define GFPOLY_HPP

#include "gf.hpp"

#include <vector>
#include <cassert>
#include <ostream>
#include <algorithm>
#include <functional>
#include <initializer_list>

namespace irrpoly {

/**
 * Класс gfpoly представляет многочлены над полем Галуа.
 * Основан на классе polynomial из библиотеки Boost.
 */
class gfpoly;

namespace {

template<typename N>
void division_impl(gfpoly &, gfpoly &, const gfpoly &, N, N);

std::pair<gfpoly, gfpoly> division(gfpoly, const gfpoly &);

std::pair<gfpoly, gfpoly> quotient_remainder(const gfpoly &, const gfpoly &);

} // namespace

class gfpoly {
private:
    gf m_field;
    std::vector<gfn> m_data;

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

    static
    std::vector<uintmax_t> into_vec(gfpoly &&poly) {
        std::vector<uintmax_t> res;
        res.reserve(poly.size());
        for (gfn &val : poly.m_data) {
            res.emplace_back(gfn::into_num(std::move(val)));
        }
        return res;
    }

private:
    gfpoly(gf field, std::vector<gfn> &&p) :
        m_field(std::move(field)), m_data(std::move(p)) {
        normalize();
    }

public:
    explicit
    gfpoly(gf field) :
        m_field(std::move(field)), m_data() {}

    gfpoly(const gf &field, const std::vector<uintmax_t> &l) :
        m_field(field), m_data() {
        m_data.reserve(l.size());
        for (uintmax_t v : l) {
            m_data.emplace_back(field, v);
        }
    }

    gfpoly(const gf &field, std::initializer_list<uintmax_t> l) :
        gfpoly(field, std::vector<uintmax_t>{l}) {}

    gfpoly(const gfpoly &p) = default;

    gfpoly(gfpoly &&p) = default;

    explicit
    gfpoly(gfn value) :
        m_field(value.field()), m_data() {
        m_data.push_back(std::move(value));
    }

    gfpoly(const gf &field, uintmax_t value) :
        m_field(field), m_data() {
        m_data.emplace_back(field, value);
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
        return m_data.size();
    }

    [[nodiscard]]
    uintmax_t degree() const {
        if (size() == 0)
            throw std::logic_error("degree() is undefined for the zero polynomial");
        return m_data.size() - 1;
    }

    gfn &operator[](uintmax_t i) {
        return m_data[i];
    }

    const gfn &operator[](uintmax_t i) const {
        return m_data[i];
    }

    [[nodiscard]]
    const std::vector<gfn> &data() const {
        return m_data;
    }

    std::vector<gfn> &data() {
        return m_data;
    }

    gfpoly &operator=(const gfpoly &p) {
        if (this != &p) {
            assert(m_field == nullptr || m_field == p.m_field);
            m_field = p.m_field;
            m_data = p.m_data;
        }
        return *this;
    }

    template<class U>
    gfpoly &operator+=(const U &value) {
        transform(value, std::plus());
        normalize();
        return *this;
    }

    template<class U>
    gfpoly &operator-=(const U &value) {
        transform(value, std::minus());
        normalize();
        return *this;
    }

    template<class U>
    gfpoly &operator*=(const U &value) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x * value; });
        normalize();
        return *this;
    }

    template<class U>
    gfpoly &operator/=(const U &value) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [&](const gfn &x) -> gfn { return x / value; });
        normalize();
        return *this;
    }

    template<class U>
    gfpoly &operator%=(const U & /*value_type*/) {
        // We can always divide by a scalar, so there is no remainder:
        this->set_zero();
        return *this;
    }

    gfpoly &operator+=(const gfpoly &value) {
        assert(field() == value.field());
        transform(value, std::plus());
        normalize();
        return *this;
    }

    gfpoly &operator-=(const gfpoly &value) {
        assert(field() == value.field());
        transform(value, std::minus());
        normalize();
        return *this;
    }

    void multiply(const gfpoly &a, const gfpoly &b) {
        assert(a.field() == b.field());
        if (!a || !b) {
            this->set_zero();
            return;
        }
        std::vector<gfn> prod(a.size() + b.size() - 1, gfn(m_field));
        for (typename std::vector<gfn>::size_type i = 0; i < a.size(); ++i)
            for (typename std::vector<gfn>::size_type j = 0; j < b.size(); ++j)
                prod[i + j] += a.m_data[i] * b.m_data[j];
        m_data.swap(prod);
    }

    gfpoly &operator*=(const gfpoly &value) {
        assert(field() == value.field());
        this->multiply(*this, value);
        return *this;
    }

    gfpoly &operator/=(const gfpoly &value) {
        assert(field() == value.field());
        *this = quotient_remainder(*this, value).first;
        return *this;
    }

    gfpoly &operator%=(const gfpoly &value) {
        assert(field() == value.field());
        *this = quotient_remainder(*this, value).second;
        return *this;
    }

    template<typename U>
    gfpoly &operator>>=(U const &n) {
        assert(n <= m_data.size());
        m_data.erase(m_data.begin(), m_data.begin() + n);
        return *this;
    }

    template<typename U>
    gfpoly &operator<<=(U const &n) {
        m_data.insert(m_data.begin(), n, gfn(m_field));
        normalize();
        return *this;
    }

    [[nodiscard]]
    bool is_zero() const {
        return m_data.empty();
    }

    explicit operator bool() const {
        return !m_data.empty();
    }

    gfpoly &set_zero() {
        m_data.clear();
        return *this;
    }

    gfpoly &normalize() {
        m_data.erase(std::find_if(
            m_data.rbegin(), m_data.rend(),
            std::not_fn(std::mem_fn(&gfn::is_zero))
        ).base(), m_data.end());
        return *this;
    }

private:
    template<class U, class R>
    gfpoly &transform(const U &value, R op) {
        assert(field() == value.field());
        if (m_data.empty())
            m_data.resize(1, gfn(m_field));
        m_data[0] = op(m_data[0], value);
        return *this;
    }

    template<class R>
    gfpoly &transform(const gfpoly &value, R op) {
        assert(field() == value.field());
        if (m_data.size() < value.size())
            m_data.resize(value.size(), gfn(m_field));
        for (uintmax_t i = 0; i < value.size(); ++i)
            m_data[i] = op(m_data[i], value[i]);
        return *this;
    }
};

gfpoly operator-(gfpoly a) {
    std::transform(a.data().begin(), a.data().end(), a.data().begin(), std::negate());
    return a;
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
    assert(a.field() == b.field());
    gfpoly result(a.field());
    result.multiply(a, b);
    return result;
}

gfpoly operator/(const gfpoly &a, const gfpoly &b) {
    assert(a.field() == b.field());
    return quotient_remainder(a, b).first;
}

gfpoly operator%(const gfpoly &a, const gfpoly &b) {
    assert(a.field() == b.field());
    return quotient_remainder(a, b).second;
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
    for (unsigned i = 0; i < poly.size(); ++i) {
        if (i)
            os << ", ";
        os << poly[i];
    }
    os << " }";
    return os;
}

namespace {
template<typename N>
void division_impl(gfpoly &q, gfpoly &u, const gfpoly &v, N n, N k) {
    assert(q.field() == u.field() && u.field() == v.field());
    q[k] = u[n + k] / v[n];
    for (N j = n + k; j > k;) {
        j--;
        u[j] -= q[k] * v[j - k];
    }
}

std::pair<gfpoly, gfpoly> division(gfpoly u, const gfpoly &v) {
    assert(u.field() == v.field() && v.size() <= u.size() && v && u);
    uintmax_t const m = u.size() - 1, n = v.size() - 1;
    uintmax_t k = m - n;
    gfpoly q(u.field());
    q.data().resize(m - n + 1, gfn(q.field()));

    do {
        division_impl(q, u, v, n, k);
    } while (k-- != 0);
    u.data().resize(n, gfn(q.field())), gfn(q.field());
    u.normalize();
    return std::make_pair(q, u);
}

/** Calculates a / b and a % b, returning the pair (quotient, remainder) together
  * because the same amount of computation yields both.
  * This function is not defined for division by zero: user beware.
  */
std::pair<gfpoly, gfpoly> quotient_remainder(const gfpoly &dividend, const gfpoly &divisor) {
    assert(dividend.field() == divisor.field() && divisor);
    if (dividend.size() < divisor.size())
        return std::make_pair(gfpoly(dividend.field()), dividend);
    return division(dividend, divisor);
}
}

} // namespace irrpoly

#endif //GFPOLY_HPP
