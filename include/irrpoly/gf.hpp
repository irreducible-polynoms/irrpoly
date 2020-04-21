/**
 * @file    gf.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include "nn.hpp"

#include <algorithm>
#include <utility>
#include <vector>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <cassert>
#include <memory>

namespace irrpoly {

/// Класс gfn представляет число в поле Галуа.
class gfn;

/// Класс gf представляет поле Галуа.
class gfbase final {
private:
    std::vector<uintmax_t> m_inv; /// Обратные элементы по умножению
    const uintmax_t m_base; /// Основание поля

    explicit
    gfbase(uintmax_t);

    friend dropbox::oxygen::nn_shared_ptr<gfbase> make_gf(uintmax_t);

public:
    [[nodiscard]]
    uintmax_t base() const;

    [[nodiscard]]
    uintmax_t mul_inv(uintmax_t) const;
};

using gf = dropbox::oxygen::nn_shared_ptr<gfbase>;

bool operator==(const gf &lb, const gf &rb) {
    return lb->base() == rb->base();
}

bool operator!=(const gf &lb, const gf &rb) {
    return lb->base() != rb->base();
}

/**
 * Класс gfn представляет число в поле GF[P]
 * @tparam P основание поля Галуа GF[P].
 */
class gfn final {
private:
    gf m_field;
    uintmax_t m_val;

public:
    /// Генерирует случайное число в пределах [0, P-1].
    static
    gfn random(const gf &field) {
        static std::random_device rd;
#ifdef __LP64__
        static std::mt19937_64 gen(rd());
#else
        static std::mt19937 gen(rd());
#endif
        std::uniform_int_distribution<uint_fast64_t> dis(0, field->base() - 1);
        return gfn(field, dis(gen));
    }

    static
    uintmax_t into_num(gfn &&val) {
        return val.m_val;
    }

    /// Конструктор по умолчанию, обнуляет переменную.
    explicit
    gfn(const gf &field) : m_field(field), m_val(0) {}

    gfn(const gf &field, const uintmax_t val) :
        m_field(field), m_val(val % m_field->base()) {}

    gfn(const gfn &other) = default;

    gfn(gfn &&other) = default;

    gfn &operator=(const gfn &other) {
        if (this != &other) {
            assert(m_field == nullptr || m_field == other.m_field);
            m_field = other.m_field;
            m_val = other.m_val;
        }
        return *this;
    }

    [[nodiscard]]
    uintmax_t base() const {
        return m_field->base();
    }

    gfn &operator=(const uintmax_t other) {
        m_val = other % base();
        return *this;
    }

    [[nodiscard]]
    const gf &field() const {
        return m_field;
    }

    /// Возвращает значение класса в виде целого числа.
    [[nodiscard]]
    uintmax_t data() const {
        return m_val;
    }

    [[nodiscard]]
    gfn operator+() const {
        return gfn{*this};
    }

    [[nodiscard]]
    gfn operator+(const gfn &other) const {
        assert(m_field == other.field());
        return gfn{m_field, m_val + other.m_val};
    }

    [[nodiscard]]
    gfn operator+(const uintmax_t other) const {
        return gfn{m_field, m_val + (other % m_field->base())};
    }

    friend
    gfn operator+(uintmax_t, const gfn &);

    gfn &operator+=(const gfn &other) {
        assert(m_field == other.field());
        m_val = (m_val + other.m_val) % base();
        return *this;
    }

    gfn operator+=(const uintmax_t other) {
        m_val = (m_val + (other % base())) % base();
        return *this;
    }

    gfn &operator++() {
        m_val = (m_val + 1) % base();
        return *this;
    }

    [[nodiscard]]
    gfn operator++(int) &{
        gfn tmp{*this};
        m_val = (m_val + 1) % base();
        return tmp;
    }

    [[nodiscard]]
    gfn operator-() const {
        gfn tmp{*this};
        tmp.m_val = base() - m_val;
        return tmp;
    }

    [[nodiscard]]
    gfn operator-(const gfn &other) const {
        assert(m_field == other.field());
        return gfn{m_field, base() + m_val - other.m_val};
    }

    [[nodiscard]]
    gfn operator-(const uintmax_t other) const {
        return gfn{m_field, base() + m_val - (other % base())};
    }

    friend
    gfn operator-(uintmax_t, const gfn &);

    gfn &operator-=(const gfn &other) {
        assert(m_field == other.field());
        m_val = (base() + m_val - other.m_val) % base();
        return *this;
    }

    gfn operator-=(const uintmax_t other) {
        m_val = (base() + m_val - (other % base())) % base();
        return *this;
    }

    gfn &operator--() {
        m_val = (base() + m_val - 1) % base();
        return *this;
    }

    [[nodiscard]]
    gfn operator--(int) &{
        gfn tmp{*this};
        m_val = (base() + m_val - 1) % base();
        return tmp;
    }

    [[nodiscard]]
    gfn operator*(const gfn &other) const {
        assert(m_field == other.field());
        return gfn{m_field, m_val * other.m_val};
    }

    [[nodiscard]]
    gfn operator*(const uintmax_t other) const {
        return gfn{m_field, m_val * (other % base())};
    }

    friend
    gfn operator*(uintmax_t, const gfn &);

    gfn &operator*=(const gfn &other) {
        assert(m_field == other.field());
        m_val = (m_val * other.m_val) % base();
        return *this;
    }

    gfn operator*=(const uintmax_t other) {
        m_val = (m_val * (other % base())) % base();
        return *this;
    }

    [[nodiscard]]
    gfn mul_inv() {
        return gfn{m_field, m_field->mul_inv(m_val)};
    }

    [[nodiscard]]
    gfn operator/(const gfn &other) const {
        assert(m_field == other.field());
        switch (other.m_val) {
        case 0:throw std::invalid_argument("division by zero");
        default:return gfn{m_field, m_val * m_field->mul_inv(other.m_val)};
        }
    }

    [[nodiscard]]
    gfn operator/(const uintmax_t other) const {
        switch (other % base()) {
        case 0:throw std::invalid_argument("division by zero");
        default:
            return gfn{m_field,
                       m_val * m_field->mul_inv(other % base())};
        }
    }

    friend
    gfn operator/(uintmax_t, const gfn &);

    gfn &operator/=(const gfn &other) {
        assert(m_field == other.field());
        switch (other.m_val) {
        case 0:throw std::invalid_argument("division by zero");
        default:m_val = (m_val * m_field->mul_inv(other.m_val)) % base();
            return *this;
        }
    }

    gfn operator/=(const uintmax_t other) {
        switch (other % base()) {
        case 0:throw std::invalid_argument("division by zero");
        default:m_val = (m_val * m_field->mul_inv(other % base())) % base();
            return *this;
        }
    }

    /// Проверяет равенство данного числа нулю.
    [[nodiscard]]
    bool is_zero() const {
        return 0 == m_val;
    }

    explicit operator bool() const {
        return !m_val;
    }

    gfn &set_zero() {
        m_val = 0;
        return *this;
    }

    bool operator==(const gfn &other) const {
        assert(m_field == other.field());
        return m_val == other.m_val;
    }

    bool operator==(const uintmax_t other) const {
        return m_val == (other % base());
    }

    bool operator!=(const gfn &other) const {
        assert(m_field == other.field());
        return m_val != other.m_val;
    }

    bool operator!=(const uintmax_t other) const {
        return m_val != (other % base());
    }

    bool operator>(const gfn &other) const {
        assert(m_field == other.field());
        return m_val > other.m_val;
    }

    bool operator>(const uintmax_t other) const {
        return m_val > (other % base());
    }

    bool operator>=(const gfn &other) const {
        assert(m_field == other.field());
        return m_val >= other.m_val;
    }

    bool operator>=(const uintmax_t other) const {
        return m_val >= (other % base());
    }

    bool operator<(const gfn &other) const {
        assert(m_field == other.field());
        return m_val < other.m_val;
    }

    bool operator<(const uintmax_t other) const {
        return m_val < (other % base());
    }

    bool operator<=(const gfn &other) const {
        assert(m_field == other.field());
        return m_val <= other.m_val;
    }

    bool operator<=(const uintmax_t other) const {
        return m_val <= (other % base());
    }

    friend
    std::ostream &operator<<(std::ostream &, const gfn &);

    friend
    std::istream &operator>>(std::istream &, gfn &);
};

[[nodiscard]]
gfn operator+(const uintmax_t other, const gfn &curr) {
    return gfn{curr.m_field, (other % curr.base()) + curr.m_val};
}

[[nodiscard]]
gfn operator-(const uintmax_t other, const gfn &curr) {
    return gfn{curr.m_field, curr.base() + (other % curr.base()) - curr.m_val};
}

[[nodiscard]]
gfn operator*(const uintmax_t other, const gfn &curr) {
    return gfn{curr.m_field, (other % curr.base()) * curr.m_val};
}

[[nodiscard]]
gfn operator/(const uintmax_t other, const gfn &curr) {
    switch (curr.m_val) {
    case 0:throw std::invalid_argument("division by zero");
    default:
        return gfn{curr.m_field,
                   (other % curr.base()) * curr.m_field->mul_inv(curr.m_val)};
    }
}

std::ostream &operator<<(std::ostream &os, const gfn &val) {
    return os << val.m_val;
}

std::istream &operator>>(std::istream &is, gfn &val) {
    is >> val.m_val;
    val.m_val %= val.base();
    return is;
}

inline
gfbase::gfbase(const uintmax_t base) : m_base(base), m_inv(base, 0) {
    if (base == 0)
        throw std::logic_error("empty field");
    if (base == 1)
        throw std::logic_error("field could contain only zero");
    if (UINTMAX_MAX / (base - 1) < (base - 1))
        throw std::logic_error("too large field");

    auto i_base = (intmax_t) base;
    auto inv_calc = [](const intmax_t base, const intmax_t val) -> uintmax_t {
        intmax_t u0 = base, u1 = 1, u2 = 0,
            v0 = val, v1 = 0, v2 = 1, w0, w1, w2, q;
        while (v0 > 0) {
            q = u0 / v0;
            w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
            u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
        }
        if (u0 > 1)
            throw std::logic_error("multiplicative inverse don't exist");
        return (uintmax_t) (u2 < 0 ? (base + u2) : (u2));
    };

    m_inv[1] = 1;
    for (uintmax_t i = 2; i < m_base; ++i) {
        if (m_inv[i])
            continue;
        m_inv[i] = inv_calc(i_base, i);
        m_inv[m_inv[i]] = i;
    }
}

[[nodiscard]]
inline
uintmax_t gfbase::base() const {
    return m_base;
}

[[nodiscard]]
inline
uintmax_t gfbase::mul_inv(const uintmax_t val) const {
    switch (val % m_base) {
    case 0:throw std::logic_error("multiplicative inverse don't exist");
    default:return m_inv[val % m_base];
    }
}

[[nodiscard]]
inline
gf make_gf(const uintmax_t base) {
    return dropbox::oxygen::nn_shared_ptr<gfbase>(dropbox::oxygen::nn(
        dropbox::oxygen::i_promise_i_checked_for_null_t{}, new gfbase(base)));
}

} // namespace irrpoly
