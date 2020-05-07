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
#include <ostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <memory>

namespace irrpoly {

/**
 * Бинарные операции над gfn корректно определены лишь в том случае, когда
 * оба числа принадлежат одному и тому же полю. По умолчанию, проверка этого
 * факта выполняется только для конфигурации Debug, чтобы в конфигурации Release
 * добиться наибольшей производительности. Однако, если в конфигурации Release
 * проверки также требуется выполнять, достаточно объявить
 * #define IRRPOLY_RELEASE_CHECKED
 * перед подключением заголовочного файла
 * #include <irrpoly.h>
 */
#if !defined(NDEBUG) || defined(IRRPOLY_RELEASE_CHECKED) // Debug or Release Checked
#define CHECK_FIELD(comparison) \
    if (!(comparison)) { \
        throw std::logic_error("field check failed"); \
    }
#else // Release
#define CHECK_FIELD(comparison)
#endif

/// Класс, содержащий основание поля Галуа и все обратные по умножению элементы.
class gfbase;

/// Поле всегда должно существовать, поэтому используется обёртка, не позволяющая
/// хранящим поле переменным иметь значение nullptr. Для него корректно определена
/// операция копирования за счёт использования примитива shared_ptr.
using gf = dropbox::oxygen::nn_shared_ptr<gfbase>;

class gfbase final {
private:
    const uintmax_t m_base; /// Основание поля
    std::vector<uintmax_t> m_inv; /// Обратные элементы по умножению

    explicit
    gfbase(uintmax_t /*base*/);

    friend auto make_gf(uintmax_t /*base*/) -> gf;

public:
    /// Возвращает основание поля Галуа.
    [[nodiscard]]
    auto base() const -> uintmax_t;

    /// Возвращает обратное по умнажению в виде числа.
    [[nodiscard]]
    auto mul_inv(uintmax_t /*val*/) const -> uintmax_t;
};

auto operator==(const gf &lb, const gf &rb) -> bool {
    return lb->base() == rb->base();
}

auto operator!=(const gf &lb, const gf &rb) -> bool {
    return lb->base() != rb->base();
}

/// Класс, представляющий число в поле Галуа с корректно заданными
/// арифметическими операциями и операторами сравнения.
class gfn final {
private:
    gf m_field; /// Поле GF[P]
    uintmax_t m_val; /// Число в поле, всегда лежит в диапазоне [0, P-1]

public:
    /// Генерирует случайное число в диапазоне [0, P-1].
    static
    auto random(const gf &field) -> gfn {
        static std::random_device rd;
#ifdef __LP64__
        static std::mt19937_64 gen(rd());
#else
        static std::mt19937 gen(rd());
#endif
        std::uniform_int_distribution<uint_fast64_t> dis(0, field->base() - 1);
        return gfn(field, dis(gen));
    }

    /// Возвращает значение класса в виде целого числа.
    [[nodiscard]]
    auto value() const -> uintmax_t {
        return m_val;
    }

    /// Конструктор по умолчанию, обнуляет переменную.
    explicit
    gfn(const gf &field) : m_field(field), m_val(0) {}

    /// Поле всегда принимается по указателю и копируется в последний момент.
    gfn(const gf &field, const uintmax_t val) :
        m_field(field), m_val(val % m_field->base()) {}

    gfn(const gfn &other) = default;

    gfn(gfn &&other) = default;

    /// Операция присваивания корректно определена только для числел из одного поля,
    /// операция перемещения также корректно определена для неинициализированных
    /// переменных (случай field = nullptr).
    auto operator=(const gfn &other) -> gfn & {
        if (this != &other) {
            CHECK_FIELD(m_field == nullptr || m_field == other.m_field)
            m_field = other.m_field;
            m_val = other.m_val;
        }
        return *this;
    }

    /// Возвращает основание поля Галуа.
    [[nodiscard]]
    auto base() const -> uintmax_t {
        return m_field->base();
    }

    auto operator=(const uintmax_t other) -> gfn & {
        m_val = other % base();
        return *this;
    }

    [[nodiscard]]
    auto field() const -> const gf & {
        return m_field;
    }

    [[nodiscard]]
    auto operator+() const -> gfn {
        return gfn{*this};
    }

    [[nodiscard]]
    auto operator+(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        return gfn{m_field, m_val + other.m_val};
    }

    [[nodiscard]]
    auto operator+(const uintmax_t other) const -> gfn {
        return gfn{m_field, m_val + (other % m_field->base())};
    }

    friend
    auto operator+(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator+=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        m_val = (m_val + other.m_val) % base();
        return *this;
    }

    auto operator+=(const uintmax_t other) -> gfn {
        m_val = (m_val + (other % base())) % base();
        return *this;
    }

    auto operator++() -> gfn & {
        m_val = (m_val + 1) % base();
        return *this;
    }

    [[nodiscard]]
    auto operator++(int) & -> gfn {
        gfn tmp{*this};
        m_val = (m_val + 1) % base();
        return tmp;
    }

    [[nodiscard]]
    auto operator-() const -> gfn {
        gfn tmp{*this};
        tmp.m_val = (base() - m_val) % base();
        return tmp;
    }

    [[nodiscard]]
    auto operator-(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        return gfn{m_field, base() + m_val - other.m_val};
    }

    [[nodiscard]]
    auto operator-(const uintmax_t other) const -> gfn {
        return gfn{m_field, base() + m_val - (other % base())};
    }

    friend
    auto operator-(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator-=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        m_val = (base() + m_val - other.m_val) % base();
        return *this;
    }

    auto operator-=(const uintmax_t other) -> gfn {
        m_val = (base() + m_val - (other % base())) % base();
        return *this;
    }

    auto operator--() -> gfn & {
        m_val = (base() + m_val - 1) % base();
        return *this;
    }

    [[nodiscard]]
    auto operator--(int) & -> gfn {
        gfn tmp{*this};
        m_val = (base() + m_val - 1) % base();
        return tmp;
    }

    [[nodiscard]]
    auto operator*(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        return gfn{m_field, m_val * other.m_val};
    }

    [[nodiscard]]
    auto operator*(const uintmax_t other) const -> gfn {
        return gfn{m_field, m_val * (other % base())};
    }

    friend
    auto operator*(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator*=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        m_val = (m_val * other.m_val) % base();
        return *this;
    }

    auto operator*=(const uintmax_t other) -> gfn {
        m_val = (m_val * (other % base())) % base();
        return *this;
    }

    /// Возвращает обратное по умножению в виде gfn.
    [[maybe_unused]] [[nodiscard]]
    auto mul_inv() -> gfn {
        return gfn{m_field, m_field->mul_inv(m_val)};
    }

    /// Операция деления определена как умножение на обратный элемент.
    [[nodiscard]]
    auto operator/(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        switch (other.m_val) {
        case 0:throw std::invalid_argument("division by zero");
        default:return gfn{m_field, m_val * m_field->mul_inv(other.m_val)};
        }
    }

    [[nodiscard]]
    auto operator/(const uintmax_t other) const -> gfn {
        switch (other % base()) {
        case 0:throw std::invalid_argument("division by zero");
        default:
            return gfn{m_field,
                       m_val * m_field->mul_inv(other % base())};
        }
    }

    friend
    auto operator/(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator/=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        switch (other.m_val) {
        case 0:throw std::invalid_argument("division by zero");
        default:m_val = (m_val * m_field->mul_inv(other.m_val)) % base();
            return *this;
        }
    }

    auto operator/=(const uintmax_t other) -> gfn {
        switch (other % base()) {
        case 0:throw std::invalid_argument("division by zero");
        default:m_val = (m_val * m_field->mul_inv(other % base())) % base();
            return *this;
        }
    }

    [[nodiscard]]
    auto is_zero() const -> bool {
        return 0 == m_val;
    }

    explicit operator bool() const {
        return 0 != m_val;
    }

    template<class charT, class traits>
    friend
    auto operator<<(std::basic_ostream<charT, traits> & /*os*/, const gfn & /*val*/)
    -> std::basic_ostream<charT, traits> &;

    template<class charT, class traits>
    friend
    auto operator>>(std::basic_istream<charT, traits> & /*is*/, gfn & /*val*/)
    -> std::basic_istream<charT, traits> &;
};

#define GFN_COMPARISON_OPERATORS(op) \
    inline \
    auto operator op(const gfn &l, const gfn &r) -> bool { \
        CHECK_FIELD(l.field() == r.field()) \
        return l.value() op r.value(); \
    } \
    \
    inline \
    auto operator op(const gfn &l, const uintmax_t r) -> bool { \
        return l.value() op (r % l.base()); \
    } \
    \
    inline \
    auto operator op(const uintmax_t l, const gfn &r) -> bool { \
        return (l % r.base()) op r.value(); \
    }

GFN_COMPARISON_OPERATORS(==)
GFN_COMPARISON_OPERATORS(!=)
GFN_COMPARISON_OPERATORS(<)
GFN_COMPARISON_OPERATORS(<=)
GFN_COMPARISON_OPERATORS(>)
GFN_COMPARISON_OPERATORS(>=)

#undef GFN_COMPARISON_OPERATORS

[[nodiscard]]
auto operator+(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, (other % curr.base()) + curr.m_val};
}

[[nodiscard]]
auto operator-(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, curr.base() + (other % curr.base()) - curr.m_val};
}

[[nodiscard]]
auto operator*(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, (other % curr.base()) * curr.m_val};
}

[[nodiscard]]
auto operator/(const uintmax_t other, const gfn &curr) -> gfn {
    switch (curr.m_val) {
    case 0:throw std::invalid_argument("division by zero");
    default:
        return gfn{curr.m_field,
                   (other % curr.base()) * curr.m_field->mul_inv(curr.m_val)};
    }
}

template<class charT, class traits>
auto operator<<(std::basic_ostream<charT, traits> &os, const gfn &val)
-> std::basic_ostream<charT, traits> & {
    return os << val.m_val;
}

template<class charT, class traits>
auto operator>>(std::basic_istream<charT, traits> &is, gfn &val)
-> std::basic_istream<charT, traits> & {
    is >> val.m_val;
    val.m_val %= val.base();
    return is;
}

/**
 * Поле обязательно должно содержать как минимум 0 и 1, кроме того, если для
 * какого-то из элементов предполагаемого поля не существует обратного элемента
 * по умножению, значит данное кольцо не является полем. При этом все операции
 * должны быть корректно определены, поэтому слишком большое поле создать
 * невозможно (в нём умножение двух максимальных элементов могло бы вызывать
 * переполнение, а значит приводить к undefined behaviour). Эта же проверка
 * гарантирует, что безнаковое значение поля может быть без дополнительных
 * проверок приведено к знаковому типу (intmax_t).
 */
inline
gfbase::gfbase(const uintmax_t base) : m_base(base), m_inv(base, 0) {
    if (base == 0) {
        throw std::logic_error("empty field");
    }
    if (base == 1) {
        throw std::logic_error("field could contain only zero");
    }
    if (UINTMAX_MAX / (base - 1) < (base - 1)) {
        throw std::logic_error("too large field");
    }

    auto i_base = static_cast<intmax_t>(base);
    auto inv_calc = [](const intmax_t base, const intmax_t val) -> uintmax_t {
        intmax_t u0 = base, u1 = 1, u2 = 0,
            v0 = val, v1 = 0, v2 = 1, w0 = 0, w1 = 0, w2 = 0, q = 0;
        while (v0 > 0) {
            q = u0 / v0;
            w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
            u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
        }
        if (u0 > 1) {
            throw std::logic_error("multiplicative inverse don't exist");
        }
        return static_cast<uintmax_t>(u2 < 0 ? (base + u2) : (u2));
    };

    m_inv[1] = 1;
    for (uintmax_t i = 2; i < m_base; ++i) {
        if (m_inv[i]) {
            continue;
        }
        m_inv[i] = inv_calc(i_base, i);
        m_inv[m_inv[i]] = i;
    }
}

[[nodiscard]]
inline
auto gfbase::base() const -> uintmax_t {
    return m_base;
}

[[nodiscard]]
inline
auto gfbase::mul_inv(const uintmax_t val) const -> uintmax_t {
    switch (val % m_base) {
    case 0:throw std::logic_error("multiplicative inverse don't exist");
    default:return m_inv[val % m_base];
    }
}

[[nodiscard]]
inline
auto make_gf(const uintmax_t base) -> gf {
    return dropbox::oxygen::nn<std::shared_ptr<gfbase>>(dropbox::oxygen::nn(
        dropbox::oxygen::i_promise_i_checked_for_null_t{}, new gfbase(base)));
}

#undef CHECK_FIELD

} // namespace irrpoly
