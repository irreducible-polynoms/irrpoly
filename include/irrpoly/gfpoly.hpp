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
#include <initializer_list>

namespace irrpoly {

/**
 * Класс gfpoly представляет многочлены над полем Галуа.
 * Основан на классе polynomial из библиотеки Boost.
 */
class gfpoly;

namespace detail {

template<typename N>
void division_impl(gfpoly &, gfpoly &, const gfpoly &, N, N);

std::pair<gfpoly, gfpoly> division(gfpoly, const gfpoly &);

/** Calculates a / b and a % b, returning the pair (quotient, remainder) together
  * because the same amount of computation yields both.
  * This function is not defined for division by zero: user beware.
  */
std::pair<gfpoly, gfpoly> quotient_remainder(const gfpoly &, const gfpoly &);

} // namespace detail

class gfpoly {
private:
    gf m_field;
    std::vector<gfn> m_data;

public:
    typedef typename std::vector<gfn>::value_type value_type;
    typedef typename std::vector<gfn>::size_type size_type;

    /// Генерирует случайный многочлен над полем GF[P] заданной степени.
    static
    gfpoly random(const gf &field, typename gfpoly::size_type degree) {
        std::vector<gfn> data;
        data.reserve(degree + 1);
        for (auto i = 0; i < degree; ++i) {
            data.emplace_back(gfn::random(field));
        }
        data.emplace_back(field, 1);
        while (data[0].is_zero()) {
            data[0] = gfn::random(field);
        }
        return gfpoly(field, std::move(data));
    }

    explicit
    gfpoly(gf field) :
        m_field(std::move(field)), m_data() {}

private:
    gfpoly(gf field, std::vector<gfn> &&p) :
        m_field(std::move(field)), m_data(std::move(p)) {
        normalize();
    }

public:
    gfpoly(gfpoly &&p) noexcept :
        m_field(std::move(p.m_field)), m_data(std::move(p.m_data)) {}

    gfpoly(const gfpoly &p) = default;

    gfpoly(const gf &field, std::initializer_list<uintmax_t> l) :
        m_field(field), m_data(make_gfn(field, l)) {}

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
    gf field() {
        return m_field;
    }

    [[nodiscard]]
    const gf &field() const {
        return m_field;
    }

    [[nodiscard]]
    size_type size() const {
        return m_data.size();
    }

    [[nodiscard]]
    size_type degree() const {
        if (size() == 0)
            throw std::logic_error("degree() is undefined for the zero polynomial");
        return m_data.size() - 1;
    }

    value_type &operator[](size_type i) {
        return m_data[i];
    }

    const value_type &operator[](size_type i) const {
        return m_data[i];
    }

    [[nodiscard]]
    const std::vector<gfn> &data() const {
        return m_data;
    }

    std::vector<gfn> &data() {
        return m_data;
    }

    gfpoly &operator=(gfpoly &&p) noexcept {
        m_field = std::move(p.m_field);
        m_data = std::move(p.m_data);
        return *this;
    }

    gfpoly &operator=(const gfpoly &p) {
        if (this != &p) {
            m_field = p.m_field;
            m_data = p.m_data;
        }
        return *this;
    }

    template<class U>
    gfpoly &operator+=(const U &value) {
        addition(value);
        normalize();
        return *this;
    }

    template<class U>
    gfpoly &operator-=(const U &value) {
        subtraction(value);
        normalize();
        return *this;
    }

    template<class U>
    gfpoly &operator*=(const U &value) {
        multiplication(value);
        normalize();
        return *this;
    }

    template<class U>
    gfpoly &operator/=(const U &value) {
        division(value);
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
        addition(value);
        normalize();
        return *this;
    }

    gfpoly &operator-=(const gfpoly &value) {
        assert(field() == value.field());
        subtraction(value);
        normalize();
        return *this;
    }

    void multiply(const gfpoly &a, const gfpoly &b) {
        assert(a.field() == b.field());
        if (!a || !b) {
            this->set_zero();
            return;
        }
        std::vector<gfn> prod(a.size() + b.size() - 1, gfn(m_field, 0));
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
        *this = detail::quotient_remainder(*this, value).first;
        return *this;
    }

    gfpoly &operator%=(const gfpoly &value) {
        assert(field() == value.field());
        *this = detail::quotient_remainder(*this, value).second;
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
        m_data.insert(m_data.begin(), n, gfn(m_field, 0));
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

    void set_zero() {
        m_data.clear();
    }

    void normalize() {
        m_data.erase(std::find_if(m_data.rbegin(), m_data.rend(),
                                  [](const gfn &x) -> bool { return !x.is_zero(); }).base(), m_data.end());
    }

private:
    template<class U, class R>
    gfpoly &addition(const U &value, R op) {
        assert(field() == value.field());
        if (m_data.empty())
            m_data.resize(1, gfn(m_field, 0));
        m_data[0] = op(m_data[0], value);
        return *this;
    }

    template<class U>
    gfpoly &addition(const U &value) {
        assert(field() == value.field());
        return addition(value, [](const auto &x, const auto &y) { return x + y; });
    }

    template<class U>
    gfpoly &subtraction(const U &value) {
        assert(field() == value.field());
        return addition(value, [](const auto &x, const auto &y) { return x - y; });
    }

    template<class R>
    gfpoly &addition(const gfpoly &value, R op) {
        assert(field() == value.field());
        if (m_data.size() < value.size())
            m_data.resize(value.size(), gfn(m_field, 0));
        for (size_type i = 0; i < value.size(); ++i)
            m_data[i] = op(m_data[i], value[i]);
        return *this;
    }

    gfpoly &addition(const gfpoly &value) {
        assert(field() == value.field());
        return addition(value, [](const auto &x, const auto &y) { return x + y; });
    }

    gfpoly &subtraction(const gfpoly &value) {
        assert(field() == value.field());
        return addition(value, [](const auto &x, const auto &y) { return x - y; });
    }

    template<class U>
    gfpoly &multiplication(const U &value) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](const gfn &x) -> gfn { return x * value; });
        return *this;
    }

    template<class U>
    gfpoly &division(const U &value) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](const gfn &x) -> gfn { return x / value; });
        return *this;
    }
};

gfpoly operator-(gfpoly a) {
    std::transform(a.data().begin(), a.data().end(), a.data().begin(), [](const auto &x) { return -x; });
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
    return detail::quotient_remainder(a, b).first;
}

gfpoly operator%(const gfpoly &a, const gfpoly &b) {
    assert(a.field() == b.field());
    return detail::quotient_remainder(a, b).second;
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

namespace detail {
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
    assert(u.field() == v.field());
    assert(v.size() <= u.size());
    assert(v);
    assert(u);

    typedef typename gfpoly::size_type N;

    N const m = u.size() - 1, n = v.size() - 1;
    N k = m - n;
    gfpoly q(u.field());
    q.data().resize(m - n + 1, gfn(q.field(), 0));

    do {
        division_impl(q, u, v, n, k);
    } while (k-- != 0);
    u.data().resize(n, gfn(q.field(), 0)), gfn(q.field(), 0);
    u.normalize();
    return std::make_pair(q, u);
}

std::pair<gfpoly, gfpoly> quotient_remainder(const gfpoly &dividend, const gfpoly &divisor) {
    assert(dividend.field() == divisor.field());
    assert(divisor);
    if (dividend.size() < divisor.size())
        return std::make_pair(gfpoly(dividend.field()), dividend);
    return detail::division(dividend, divisor);
}
}

} // namespace irrpoly

#endif //GFPOLY_HPP
