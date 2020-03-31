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
#include <utility>
#include <vector>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <cassert>
#include <memory>

#include <array>

namespace irrpoly {

    //Класс gfn представляет число в поле Галуа.
    class gfn;

    /// Класс gf представляет поле Галуа.
    class gf final {
    public:
        using value_type = uintmax_t;

    private:
        std::vector<value_type> m_inv; /// Обратные элементы по умножению
        const value_type m_base; /// Основание поля

        explicit
        gf(uintmax_t);

        friend std::shared_ptr<gf> make_gf(uintmax_t);

    public:
        [[nodiscard]]
        value_type base() const;

        [[nodiscard]]
        value_type mul_inv(value_type) const;
    };

    using gfp = std::shared_ptr<gf>;

    namespace detail {

        /**
         * @param base основание поля Галуа GF[base]
         * @param val элемент поля GF[base], для которого требуется найти обратный
         * @return обратный элемент по умножению в поле GF[base]
         */
        [[nodiscard]]
        typename gf::value_type inv_calc(const uintmax_t base, const typename gf::value_type val) {
            intmax_t u0 = static_cast<intmax_t>(base), u1 = 1, u2 = 0,
                     v0 = val, v1 = 0, v2 = 1, w0 = 0, w1 = 0, w2 = 0, q = 0;
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
                    if (u0 > 1) throw std::logic_error("multiplicative inverse don't exist");
                    return u2 < 0 ?
                        static_cast<typename gf::value_type>(base + u2) :
                        static_cast<typename gf::value_type>(u2);
            }
        }

    }

    /**
     * Класс gfn представляет число в поле GF[P]
     * @tparam P основание поля Галуа GF[P].
     */
    class gfn final {
    public:
        using value_type = gf::value_type;

    private:
        gfp m_field;
        value_type m_val;

    public:
        /// Генерирует случайное число в пределах [0, P-1].
        static
        gfn random(gfp field) {
            static std::random_device rd;
#ifdef __LP64__
            static std::mt19937_64 gen(rd());
#else
            static std::mt19937 gen(rd());
#endif
            std::uniform_int_distribution<uint_fast64_t> dis(0,  field->base() - 1);
            return gfn(std::move(field), dis(gen));
        }

        /// Конструктор по умолчанию, обнуляет переменную.
        explicit
        gfn(gfp field) : m_field(std::move(field)), m_val(0) {}

        gfn(gfp field, const value_type val) : m_field(std::move(field)), m_val(val % m_field->base()) {}

        gfn(const gfn &other) = default;

        gfn(gfn &&other) = default;

        gfn &operator=(const gfn &other) {
            if (this != &other) {
                gfn tmp {other.m_field, other.m_val};
                std::swap(*this, tmp);
            }
            return *this;
        }

        gfn &operator=(gfn &&other) noexcept {
            m_field = std::move(other.m_field);
            m_val = other.m_val;
            return *this;
        }

        [[nodiscard]]
        const gfp &field() const {
            return m_field;
        }

        [[nodiscard]]
        gfp field() {
            return m_field;
        }

        /// Возвращает значение класса в виде целого числа.
        [[nodiscard]]
        value_type data() const {
            return m_val;
        }

        [[nodiscard]]
        gfn operator+() const {
            return gfn {*this};
        }

        [[nodiscard]]
        gfn operator+(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return gfn {m_field, m_val + other.m_val};
        }

        [[nodiscard]]
        gfn operator+(const value_type other) const {
            return gfn {m_field, m_val + (other % m_field->base())};
        }

        friend
        gfn operator+(value_type, const gfn &);

        gfn &operator+=(const gfn &other) {
            assert(m_field->base() == other.m_field->base());
            m_val = (m_val + other.m_val) % m_field->base();
            return *this;
        }

        [[nodiscard]]
        gfn operator+=(const value_type other) {
            m_val = (m_val + (other % m_field->base())) % m_field->base();
            return *this;
        }

        gfn &operator++() {
            m_val = (m_val + 1) % m_field->base();
            return *this;
        }

        [[nodiscard]]
        gfn operator++(int) & {
            gfn tmp {*this};
            m_val = (m_val + 1) % m_field->base();
            return tmp;
        }

        [[nodiscard]]
        gfn operator-() const {
            gfn tmp {*this};
            tmp.m_val = m_field->base() - m_val;
            return tmp;
        }

        [[nodiscard]]
        gfn operator-(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return gfn {m_field, m_field->base() + m_val - other.m_val};
        }

        [[nodiscard]]
        gfn operator-(const value_type other) const {
            return gfn {m_field, m_field->base() + m_val - (other % m_field->base())};
        }

        friend
        gfn operator-(value_type, const gfn &);

        gfn &operator-=(const gfn &other) {
            assert(m_field->base() == other.m_field->base());
            m_val = (m_field->base() + m_val - other.m_val) % m_field->base();
            return *this;
        }

        [[nodiscard]]
        gfn operator-=(const value_type other) {
            m_val = (m_field->base() + m_val - (other % m_field->base())) % m_field->base();
            return *this;
        }

        gfn &operator--() {
            m_val = (m_field->base() + m_val - 1) % m_field->base();
            return *this;
        }

        [[nodiscard]]
        gfn operator--(int) & {
            gfn tmp {*this};
            m_val = (m_field->base() + m_val - 1) % m_field->base();
            return tmp;
        }

        [[nodiscard]]
        gfn operator*(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return gfn {m_field, m_val * other.m_val};
        }

        [[nodiscard]]
        gfn operator*(const value_type other) const {
            return gfn {m_field, m_val * (other % m_field->base())};
        }

        friend
        gfn operator*(value_type, const gfn &);
        
        gfn &operator*=(const gfn &other) {
            assert(m_field->base() == other.m_field->base());
            m_val = (m_val * other.m_val) % m_field->base();
            return *this;
        }

        [[nodiscard]]
        gfn operator*=(const value_type other) {
            m_val = (m_val * (other % m_field->base())) % m_field->base();
            return *this;
        }
        
        [[nodiscard]]
        gfn operator/(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            switch (other.m_val) {
                case 0:
                    throw std::invalid_argument("division by zero");
                default:
                    return gfn {m_field, m_val * m_field->mul_inv(other.m_val)};
            }
        }

        [[nodiscard]]
        gfn operator/(const value_type other) const {
            switch (other % m_field->base()) {
                case 0:
                    throw std::invalid_argument("division by zero");
                default:
                    return gfn {m_field, m_val * m_field->mul_inv(other % m_field->base())};
            }
        }

        friend
        gfn operator/(value_type, const gfn &);
        
        gfn &operator/=(const gfn &other) {
            assert(m_field->base() == other.m_field->base());
            switch (other.m_val) {
                case 0:
                    throw std::invalid_argument("division by zero");
                default:
                    m_val = (m_val * m_field->mul_inv(other.m_val)) % m_field->base();
                    return *this;
            }
        }

        [[nodiscard]]
        gfn operator/=(const value_type other) {
            switch (other % m_field->base()) {
                case 0:
                    throw std::invalid_argument("division by zero");
                default:
                    m_val = (m_val * m_field->mul_inv(other % m_field->base())) % m_field->base();
                    return *this;
            }
        }

        /// Проверяет равенство данного числа нулю.
        [[nodiscard]]
        bool is_zero() const {
            return 0 == m_val;
        }
        
        bool operator==(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return m_val == other.m_val;
        }
        
        bool operator==(const value_type other) const {
            return m_val == (other % m_field->base());
        }
        
        bool operator!=(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return m_val != other.m_val;
        }
        
        bool operator!=(const value_type other) const {
            return m_val != (other % m_field->base());
        }
        
        bool operator>(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return m_val > other.m_val;
        }
        
        bool operator>(const value_type other) const {
            return m_val > (other % m_field->base());
        }
        
        bool operator>=(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return m_val >= other.m_val;
        }
        
        bool operator>=(const value_type other) const {
            return m_val >= (other % m_field->base());
        }
        
        bool operator<(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return m_val < other.m_val;
        }
        
        bool operator<(const value_type other) const {
            return m_val < (other % m_field->base());
        }
        
        bool operator<=(const gfn &other) const {
            assert(m_field->base() == other.m_field->base());
            return m_val <= other.m_val;
        }
        
        bool operator<=(const value_type other) const {
            return m_val <= (other % m_field->base());
        }

        friend
        std::ostream &operator<<(std::ostream &, const gfn &);

        friend
        std::istream &operator>>(std::istream &, gfn &);
    };
    
    [[nodiscard]]
    gfn operator+(const typename gfn::value_type other, const gfn &curr) {
        return gfn {curr.m_field, (other % curr.m_field->base()) + curr.m_val};
    }
    
    [[nodiscard]]
    gfn operator-(const typename gfn::value_type other, const gfn &curr) {
        return gfn {curr.m_field, curr.m_field->base() + (other % curr.m_field->base()) - curr.m_val};
    }
    
    [[nodiscard]]
    gfn operator*(const typename gfn::value_type other, const gfn &curr) {
        return gfn {curr.m_field, (other % curr.m_field->base()) * curr.m_val};
    }
    
    [[nodiscard]]
    gfn operator/(const typename gfn::value_type other, const gfn &curr) {
        switch (curr.m_val) {
            case 0:
                throw std::invalid_argument("division by zero");
            default:
                return gfn {curr.m_field, (other % curr.m_field->base()) * curr.m_field->mul_inv(curr.m_val)};
        }
    }

    std::ostream &operator<<(std::ostream &os, const gfn &val) {
        return os << val.m_val;
    }

    std::istream &operator>>(std::istream &is, gfn &val) {
        is >> val.m_val;
        val.m_val %= val.m_field->base();
        return is;
    }

    inline
    gf::gf(const uintmax_t base) : m_base(base), m_inv(base, 0) {
        assert(base > 0);
        for (value_type i = 1; i < m_base; ++i) {
            if (m_inv[i]) continue;
            m_inv[i] = detail::inv_calc(m_base, i);
            m_inv[m_inv[i]] = i;
        }
    }

    [[nodiscard]]
    inline
    typename gf::value_type gf::base() const {
        return m_base;
    }

    [[nodiscard]]
    inline
    typename gf::value_type gf::mul_inv(const typename gf::value_type val) const {
        switch(val % m_base) {
            case 0:
                throw std::logic_error("multiplicative inverse don't exist");
            default:
                return m_inv[val % m_base];
        }
    }

    [[nodiscard]]
    inline
    gfp make_gf(const uintmax_t base) {
        return std::shared_ptr<gf>(new gf(base));
    }

    [[nodiscard]]
    inline
    gfn make_gfn(gfp field, const typename gf::value_type val) {
        return gfn {std::move(field), val};
    }

    [[nodiscard]]
    inline
    std::vector<gfn> make_gfn(const gfp &field, const std::vector<typename gf::value_type>& val) {
        std::vector<gfn> res;
        res.reserve(val.size());
        for (typename gf::value_type v : val) {
            res.emplace_back(field, v);
        }
        return res;
    }

    [[nodiscard]]
    inline
    std::vector<gfn> make_gfn(const gfp &field, const std::initializer_list<typename gf::value_type> val) {
        return make_gfn(field, std::vector<typename gf::value_type>{val});
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
        std::ostream &operator<<(std::ostream &, const gf_old<Q> &);

        template<uint32_t Q>
        friend
        std::istream &operator>>(std::istream &, gf_old<Q> &);
    };

    /// Генерирует случайное число в пределах [0, P-1].
    template<uint32_t P>
    gf_old<P> gf_old<P>::random() noexcept {
        static std::random_device rd;
        static std::mt19937_64 gen(rd());
        static std::uniform_int_distribution<uint_fast64_t> dis(0, P - 1);
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
        static std::array<gf_old<P>, P> arr{};
        int_fast32_t u0 = P, u1 = 1, u2 = 0, v0 = v, v1 = 0, v2 = 1, w0, w1, w2, q;
        switch (v) {
            case 1:
                return gf_old<P> {*this };
            case 0:
                throw std::logic_error("multiplicative inverse don't exist");
            default:
                switch (arr[v].v) {
                    case 1: // marked as irreversible
                        throw std::logic_error("multiplicative inverse don't exist");
                    case 0: // not calculated
                        while (v0) {
                            q = u0 / v0;
                            w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
                            u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
                        }
                        if (u0 > 1) {
                            arr[v].v = 1; // mark as irreversible
                            throw std::logic_error("multiplicative inverse don't exist");
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
        if (val.v == 0) { throw std::invalid_argument("division by zero"); }
        return gf_old<P> {v * val.mul_inv().v};
    }

    template<uint32_t P>
    gf_old<P> &gf_old<P>::operator/=(const gf_old<P> &val) noexcept(false) {
        if (val.v == 0) { throw std::invalid_argument("division by zero"); }
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
    std::ostream &operator<<(std::ostream &os, const gf_old<P> &val) {
        return os << val.v;
    }

    template<uint32_t P>
    std::istream &operator>>(std::istream &is, gf_old<P> &val) {
        return is >> val.v;
    }

} // namespace irrpoly

#endif //GF_HPP