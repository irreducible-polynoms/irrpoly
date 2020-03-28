/**
 * @file    gf.hpp
 * @author  Vadim Piven <vadim@piven.tech>, Zimin Fedor <zimfv@yandex.ru>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#ifndef GF_HPP
#define GF_HPP

#include <algorithm>
#include <array>
#include <vector>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <cassert>

namespace irrpoly {

    namespace detail {

        /**
         * @param base основание поля Галуа GF[base]
         * @param val элемент поля GF[base], для которого требуется найти обратный
         * @return обратный элемент по умножению в поле GF[base]
         */
        template<typename T>
        [[nodiscard]]
        constexpr
        T inv_calc(const T base, const T val) noexcept(false) {
            intmax_t u0 = base, u1 = 1, u2 = 0, v0 = val, v1 = 0, v2 = 1, w0 = 0, w1 = 0, w2 = 0, q = 0;
            switch (val) {
                case 1:
                    return 1;
                case 0:
                    throw std::logic_error("multiplicative inverse don't exist");
                default:
                    while (v0 > 0) {
                        q = u0 / v0;
                        w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
                        u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
                    }
                    if (u0 > 1) throw ::std::logic_error("multiplicative inverse don't exist");
                    return u2 < 0 ? base + u2 : u2;
            }
        }

    }

    template <uintmax_t P>
    class gfn;

    /**
     * Класс gf представляет из себя поле GF[P].
     * @tparam P основание поля Галуа GF[P].
     */
    template<uintmax_t P = 0>
    class gf final {
    public:
        using gf_type = typename std::enable_if<UINTMAX_MAX / (P - 1) >= (P - 1),
            typename std::conditional<UINT_FAST8_MAX / (P - 1) >= (P - 1), uint_fast8_t,
                    typename std::conditional<UINT_FAST16_MAX / (P - 1) >= (P - 1), uint_fast16_t,
                            typename std::conditional<UINT_FAST32_MAX / (P - 1) >= (P - 1), uint_fast32_t,
                                    uintmax_t>::type
                    >::type
            >::type
        >::type;

    private:
        std::array<gf_type, P> m_inv; /// Обратные элементы по умножению

        template<gf_type... k>
        constexpr auto make_inv_helper(std::integer_sequence<gf_type, k...>)
        -> std::array<gf_type, sizeof...(k)> {
            return {{ [&]{
                if constexpr (k > 0)
                    return detail::inv_calc<gf_type>(sizeof...(k), k);
                else return 0;
            }()... }};
        }

    public:
        explicit constexpr
        gf(const uintmax_t base = P) noexcept(false)
            : m_inv(make_inv_helper(std::make_integer_sequence<gf_type, P>())) {
            assert(base == P);
        }

        [[nodiscard]]
        constexpr
        gf_type base() const noexcept {
            return P;
        }

        [[nodiscard]]
        constexpr
        gf_type mul_inv(const gf_type val) const noexcept(false) {
            assert(val < P);
            switch(val) {
                case 0:
                    throw std::logic_error("multiplicative inverse don't exist");
                default:
                    return m_inv[val];
            }
        }
    };

    template<>
    class gf<0> final {
    public:
        using gf_type = uintmax_t;

    private:
        std::vector<gf_type> m_inv; /// Обратные элементы по умножению
        const gf_type m_base; /// Основание поля

    public:
        explicit
        gf(const uintmax_t base) noexcept(false) : m_base(base), m_inv(base, 0) {
            assert(m_base > 0);
            for (gf_type i = 1; i < m_base; ++i) {
                if (m_inv[i]) continue;
                m_inv[i] = detail::inv_calc<gf_type>(m_base, i);
                m_inv[m_inv[i]] = i;
            }
        }

        [[nodiscard]]
        constexpr
        gf_type base() const noexcept {
            return m_base;
        }

        [[nodiscard]]
        constexpr
        gf_type mul_inv(const gf_type val) const noexcept(false) {
            assert(val < m_base);
            switch(val) {
                case 0:
                    throw std::logic_error("multiplicative inverse don't exist");
                default:
                    return m_inv[val];
            }
        }
    };

    template <uintmax_t P = 0>
    class gfn final {
    public:
        using gf_type = typename gf<P>::gf_type;

    private:
        const gf<P> &m_field;
        gf_type m_val;

    public:
        /// Генерирует случайное число в пределах [0, P-1].
        gfn random() noexcept {
            static std::random_device rd;
#ifdef __LP64__
            static std::mt19937_64 gen(rd());
#else
            static std::mt19937 gen(rd());
#endif
            std::uniform_int_distribution<uint_fast64_t> dis(0,  m_field.base() - 1);
            return gfn<P>(dis(gen));
        }

        /// Конструктор по умолчанию, обнуляет переменную.
        explicit constexpr
        gfn(const gf<P> &field) noexcept : m_field(field), m_val(0) {}

        constexpr
        gfn(const gf<P> &field, const gf_type val) noexcept : m_field(field), m_val(val % field.base()) {}

        template<uintmax_t Q>
        explicit constexpr
        gfn(const gfn<Q> &other) noexcept : m_field(other.m_field), m_val(other.m_val) {
            static_assert(0 == P || Q == P);
        }

        template<uintmax_t Q>
        constexpr
        gfn &operator=(const gfn<Q> &other) noexcept {
            static_assert(0 == P || Q == P);
            if (this != &other) {
                m_field = other.m_field;
                m_val = other.m_val;
            }
            return *this;
        }

        [[nodiscard]]
        constexpr
        const gf<P>& field() const noexcept {
            return m_field;
        }

        /// Возвращает значение класса в виде целого числа.
        [[nodiscard]]
        constexpr
        gf_type data() const noexcept {
            return m_val;
        }

        [[nodiscard]]
        constexpr
        gfn operator+() const noexcept {
            return gfn {*this};
        }

        template<uintmax_t Q>
        [[nodiscard]]
        constexpr
        gfn operator+(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return gfn {m_field, m_val + other.m_val};
        }

        template<uintmax_t Q>
        constexpr
        gfn &operator+=(const gfn<Q> &other) noexcept {
            assert(m_field.base() == other.m_field.base());
            m_val = (m_val + other.m_val) % m_field.base();
            return *this;
        }

        constexpr
        gfn &operator++() noexcept {
            m_val = (m_val + 1) % m_field.base();
            return *this;
        }

        [[nodiscard]]
        constexpr
        gfn operator++(int) & noexcept {
            gfn tmp {*this};
            m_val = (m_val + 1) % m_field.base();
            return tmp;
        }

        [[nodiscard]]
        constexpr
        gfn operator-() const noexcept {
            gfn tmp {*this};
            tmp.m_val = m_field.base() - m_val;
            return tmp;
        }

        template<uintmax_t Q>
        [[nodiscard]]
        constexpr
        gfn operator-(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return gfn {m_field, m_field.base() + m_val - other.m_val};
        }

        template<uintmax_t Q>
        constexpr
        gfn &operator-=(const gfn<Q> &other) noexcept {
            assert(m_field.base() == other.m_field.base());
            m_val = (m_field.base() + m_val - other.m_val) % m_field.base();
            return *this;
        }

        constexpr
        gfn &operator--() noexcept {
            m_val = (m_field.base() + m_val - 1) % m_field.base();
            return *this;
        }

        [[nodiscard]]
        constexpr
        gfn operator--(int) & noexcept {
            gfn tmp {*this};
            m_val = (m_field.base() + m_val - 1) % m_field.base();
            return tmp;
        }

        template<uintmax_t Q>
        [[nodiscard]]
        constexpr
        gfn operator*(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return gfn {m_field, m_val * other.m_val};
        }

        template<uintmax_t Q>
        constexpr
        gfn &operator*=(const gfn<Q> &other) noexcept {
            assert(m_field.base() == other.m_field.base());
            m_val = (m_val * other.m_val) % m_field.base();
            return *this;
        }

        template<uintmax_t Q>
        [[nodiscard]]
        constexpr
        gfn operator/(const gfn<Q> &other) const noexcept(false) {
            assert(m_field.base() == other.m_field.base());
            if (other.m_val == 0) { throw ::std::invalid_argument("division by zero"); }
            return gfn {m_field, m_val * m_field.mul_inv(other.m_val)};
        }

        template<uintmax_t Q>
        constexpr
        gfn &operator/=(const gfn<Q> &other) noexcept(false) {
            assert(m_field.base() == other.m_field.base());
            if (other.m_val == 0) { throw ::std::invalid_argument("division by zero"); }
            m_val = (m_val * m_field.mul_inv(other.m_val)) % m_field.base();
            return *this;
        }

        /// Проверяет равенство данного числа нулю.
        [[nodiscard]]
        constexpr
        bool is_zero() const noexcept {
            return m_val == 0;
        }

        template<uintmax_t Q>
        constexpr
        bool operator==(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return m_val == other.m_val;
        }

        template<uintmax_t Q>
        constexpr
        bool operator!=(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return m_val != other.m_val;
        }

        template<uintmax_t Q>
        constexpr
        bool operator>(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return m_val > other.m_val;
        }

        template<uintmax_t Q>
        constexpr
        bool operator>=(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return m_val >= other.m_val;
        }

        template<uintmax_t Q>
        constexpr
        bool operator<(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return m_val < other.m_val;
        }

        template<uintmax_t Q>
        constexpr
        bool operator<=(const gfn<Q> &other) const noexcept {
            assert(m_field.base() == other.m_field.base());
            return m_val <= other.m_val;
        }

        friend
        std::ostream &operator<<(std::ostream &, const gfn &);

        friend
        std::istream &operator>>(std::istream &, gfn &);
    };

    template<uintmax_t P>
    std::ostream &operator<<(std::ostream &os, const gfn<P> &val) {
        return os << val.m_val;
    }

    template<uintmax_t P>
    std::istream &operator>>(std::istream &is, gfn<P> &val) {
        is >> val.m_val;
        val.m_val %= val.m_field.base();
        return is;
    }

    /****************************************************************************************/

    /**
     * Класс gf представляет из себя число над полем GF[P].
     * @tparam P основание поля Галуа GF[P], по умолчанию рассматривается поле GF[2].
     * Существование поля Галуа для заданного P не гарантируется, за использование корректного
     * значения для P несут ответственность пользователи данного класса.
     */
    template<uint32_t P = 2>
    class gf_old {
    public:
        using gf_type = uint_fast64_t;

    private:
        gf_type v;

    public:
        static
        gf_old<P> random() noexcept;

        constexpr
        gf_old() noexcept;

        constexpr
        gf_old(gf_type) noexcept;

        constexpr
        gf_old(const gf_old<P> &) noexcept;

        constexpr
        gf_old<P> &operator=(const gf_old<P> &) noexcept;

        [[nodiscard]]
        constexpr
        gf_type data() const noexcept;

        [[nodiscard]]
        constexpr
        gf_old<P> operator+() const noexcept;

        [[nodiscard]]
        constexpr
        gf_old<P> operator+(const gf_old<P> &) const noexcept;

        constexpr
        gf_old<P> &operator+=(const gf_old<P> &) noexcept;

        constexpr
        gf_old<P> &operator++() noexcept;

        [[nodiscard]]
        constexpr
        gf_old<P> operator++(int) & noexcept;

        [[nodiscard]]
        constexpr
        gf_old<P> operator-() const noexcept;

        [[nodiscard]]
        constexpr
        gf_old<P> operator-(const gf_old<P> &) const noexcept;

        constexpr
        gf_old<P> &operator-=(const gf_old<P> &) noexcept;

        constexpr
        gf_old<P> &operator--() noexcept;

        [[nodiscard]]
        constexpr
        gf_old<P> operator--(int) & noexcept;

        [[nodiscard]]
        constexpr
        gf_old<P> operator*(const gf_old<P> &) const noexcept;

        constexpr
        gf_old<P> &operator*=(const gf_old<P> &) noexcept;

        [[nodiscard]]
        gf_old<P> mul_inv() const noexcept(false);

        [[nodiscard]]
        gf_old<P> operator/(const gf_old<P> &) const noexcept(false);

        gf_old<P> &operator/=(const gf_old<P> &) noexcept(false);

        [[nodiscard]]
        constexpr
        bool is_zero() const noexcept;

        constexpr
        bool operator==(const gf_old<P> &) const noexcept;

        constexpr
        bool operator!=(const gf_old<P> &) const noexcept;

        constexpr
        bool operator>(const gf_old<P> &) const noexcept;

        constexpr
        bool operator>=(const gf_old<P> &) const noexcept;

        constexpr
        bool operator<(const gf_old<P> &) const noexcept;

        constexpr
        bool operator<=(const gf_old<P> &) const noexcept;

        template<uint32_t Q>
        friend
        ::std::ostream &operator<<(::std::ostream &, const gf_old<Q> &);

        template<uint32_t Q>
        friend
        ::std::istream &operator>>(::std::istream &, gf_old<Q> &);
    };

    /// Генерирует случайное число в пределах [0, P-1].
    template<uint32_t P>
    gf_old<P> gf_old<P>::random() noexcept {
        static ::std::random_device rd;
        static ::std::mt19937_64 gen(rd());
        static ::std::uniform_int_distribution<uint_fast64_t> dis(0, P - 1);
        return gf_old<P>(dis(gen));
    }

    /// Конструктор по умолчанию, обнуляет переменную.
    template<uint32_t P>
    constexpr
    gf_old<P>::gf_old() noexcept : v(0) {}

    template<uint32_t P>
    constexpr
    gf_old<P>::gf_old(const gf_type val) noexcept : v(val % P) {}

    template<uint32_t P>
    constexpr
    gf_old<P>::gf_old(const gf_old<P> &val) noexcept : v(val.v) {}

    template<uint32_t P>
    constexpr
    gf_old<P> &gf_old<P>::operator=(const gf_old<P> &val) noexcept {
        v = val.v;
        return *this;
    }

    /// Возвращает значение класса в виде целого числа.
    template<uint32_t P>
    [[nodiscard]]
    constexpr
    typename gf_old<P>::gf_type gf_old<P>::data() const noexcept {
        return v;
    }

    template<uint32_t P>
    [[nodiscard]]
    constexpr
    gf_old<P> gf_old<P>::operator+() const noexcept {
        return gf_old<P> {*this};
    }

    template<uint32_t P>
    [[nodiscard]]
    constexpr
    gf_old<P> gf_old<P>::operator+(const gf_old<P> &val) const noexcept {
        return gf_old<P> {v + val.v};
    }

    template<uint32_t P>
    constexpr
    gf_old<P> &gf_old<P>::operator+=(const gf_old<P> &val) noexcept {
        v = (v + val.v) % P;
        return *this;
    }

    template<uint32_t P>
    constexpr
    gf_old<P> &gf_old<P>::operator++() noexcept {
        v = (v + 1) % P;
        return *this;
    }

    template<uint32_t P>
    [[nodiscard]]
    constexpr
    gf_old<P> gf_old<P>::operator++(int) & noexcept {
        gf_old<P> tmp {*this};
        v = (v + 1) % P;
        return tmp;
    }

    template<uint32_t P>
    [[nodiscard]]
    constexpr
    gf_old<P> gf_old<P>::operator-() const noexcept {
        gf_old<P> tmp {*this};
        tmp.v = P - v;
        return tmp;
    }

    template<uint32_t P>
    [[nodiscard]]
    constexpr
    gf_old<P> gf_old<P>::operator-(const gf_old<P> &val) const noexcept {
        return gf_old<P> {(P + v - val.v) % P};
    }

    template<uint32_t P>
    constexpr
    gf_old<P> &gf_old<P>::operator-=(const gf_old<P> &val) noexcept {
        v = (P + v - val.v) % P;
        return *this;
    }

    template<uint32_t P>
    constexpr
    gf_old<P> &gf_old<P>::operator--() noexcept {
        v = (P + v - 1) % P;
        return *this;
    }

    template<uint32_t P>
    [[nodiscard]]
    constexpr
    gf_old<P> gf_old<P>::operator--(int) & noexcept {
        gf_old<P> tmp {*this };
        v = (P + v - 1) % P;
        return tmp;
    }

    template<uint32_t P>
    [[nodiscard]]
    constexpr
    gf_old<P> gf_old<P>::operator*(const gf_old<P> &val) const noexcept {
        return gf_old<P> {v * val.v};
    }

    template<uint32_t P>
    constexpr
    gf_old<P> &gf_old<P>::operator*=(const gf_old<P> &val) noexcept {
        v = (v * val.v) % P;
        return *this;
    }

    /**
     * Находит обратный по умножению (multiplicative inverse) элемент для данного элемента GF[P].
     * Найденное однажды значение заносится в массив и повторно не вычисляется.
     * Для поиска используется расширенный алгоритм Евклида, реализация с минимальными
     * модификациями копирует код из библиотеки Boost 1.71.0.
     */
    template<uint32_t P>
    [[nodiscard]]
    gf_old<P> gf_old<P>::mul_inv() const noexcept(false) {
        static ::std::array<gf_old<P>, P> arr{};
        int_fast32_t u0 = P, u1 = 1, u2 = 0, v0 = v, v1 = 0, v2 = 1, w0, w1, w2, q;
        switch (v) {
            case 1:
                return gf_old<P> {*this };
            case 0:
                throw ::std::logic_error("multiplicative inverse don't exist");
            default:
                switch (arr[v].v) {
                    case 1: // marked as irreversible
                        throw ::std::logic_error("multiplicative inverse don't exist");
                    case 0: // not calculated
                        while (v0) {
                            q = u0 / v0;
                            w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
                            u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
                        }
                        if (u0 > 1) {
                            arr[v].v = 1; // mark as irreversible
                            throw ::std::logic_error("multiplicative inverse don't exist");
                        }
                        arr[v].v = u2 < 0 ? P + u2 : u2; // calculated is inverse for this
                        arr[arr[v].v].v = v; // this is inverse for calculated
                    default:
                        return arr[v];
                }
        }
    }

    /**
     * Деление реализуется как умножение на обратный по умножению.
     * Если обратный не существует - деление не возможно.
     * Для существование обратного необходимо, чтобы GF[P] было полем (зависит от выбора P).
     */
    template<uint32_t P>
    [[nodiscard]]
    gf_old<P> gf_old<P>::operator/(const gf_old<P> &val) const noexcept(false) {
        if (val.v == 0) { throw ::std::invalid_argument("division by zero"); }
        return gf_old<P> {v * val.mul_inv().v};
    }

    template<uint32_t P>
    gf_old<P> &gf_old<P>::operator/=(const gf_old<P> &val) noexcept(false) {
        if (val.v == 0) { throw ::std::invalid_argument("division by zero"); }
        v = (v * val.mul_inv().v) % P;
        return *this;
    }

    /// Проверяет равенство данного числа нулю.
    template<uint32_t P>
    [[nodiscard]]
    constexpr
    bool gf_old<P>::is_zero() const noexcept {
        return !v;
    }

    template<uint32_t P>
    constexpr
    bool gf_old<P>::operator==(const gf_old<P> &val) const noexcept {
        return v == val.v;
    }

    template<uint32_t P>
    constexpr
    bool gf_old<P>::operator!=(const gf_old<P> &val) const noexcept {
        return v != val.v;
    }

    template<uint32_t P>
    constexpr
    bool gf_old<P>::operator>(const gf_old<P> &val) const noexcept {
        return v > val.v;
    }

    template<uint32_t P>
    constexpr
    bool gf_old<P>::operator>=(const gf_old<P> &val) const noexcept {
        return v >= val.v;
    }

    template<uint32_t P>
    constexpr
    bool gf_old<P>::operator<(const gf_old<P> &val) const noexcept {
        return v < val.v;
    }

    template<uint32_t P>
    constexpr
    bool gf_old<P>::operator<=(const gf_old<P> &val) const noexcept {
        return v <= val.v;
    }

    template<uint32_t P>
    ::std::ostream &operator<<(::std::ostream &os, const gf_old<P> &val) {
        return os << val.v;
    }

    template<uint32_t P>
    ::std::istream &operator>>(::std::istream &is, gf_old<P> &val) {
        return is >> val.v;
    }

} // namespace irrpoly

#endif //GF_HPP