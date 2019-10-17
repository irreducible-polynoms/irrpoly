#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>
#include <cassert>
#include <ostream>
#include <algorithm>
#include <initializer_list>

template<typename T>
class Polynomial;

namespace detail {
/**
* Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
* Chapter 4.6.1, Algorithm D: Division of polynomials over a field.
*
* @tparam  T   Coefficient type, must be not be an integer.
*
* Template-parameter T actually must be a field but we don't currently have that
* subtlety of distinction.
*/
    template<typename T, typename N>
    void
    division_impl(Polynomial<T> &q, Polynomial<T> &u, const Polynomial<T> &v, N n, N k) {
        q[k] = u[n + k] / v[n];
        for (N j = n + k; j > k;) {
            j--;
            u[j] -= q[k] * v[j - k];
        }
    }

    template<class T, class N>
    T integer_power(T t, N n) {
        switch (n) {
            case 0:
                return static_cast<T>(1u);
            case 1:
                return t;
            case 2:
                return t * t;
            case 3:
                return t * t * t;
        }
        T result = integer_power(t, n / 2);
        result *= result;
        if (n & 1)
            result *= t;
        return result;
    }
    
/**
* Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
* Chapter 4.6.1, Algorithm D and R: Main loop.
*
* @param   u   Dividend.
* @param   v   Divisor.
*/
    template<typename T>
    std::pair<Polynomial<T>, Polynomial<T> >
    division(Polynomial<T> u, const Polynomial<T> &v) {
        assert(v.size() <= u.size());
        assert(v);
        assert(u);

        typedef typename Polynomial<T>::size_type N;

        N const m = u.size() - 1, n = v.size() - 1;
        N k = m - n;
        Polynomial<T> q;
        q.data().resize(m - n + 1);

        do {
            division_impl(q, u, v, n, k);
        } while (k-- != 0);
        u.data().resize(n);
        u.normalize(); // Occasionally, the remainder is zeroes.
        return std::make_pair(q, u);
    }

//
// These structures are the same as the void specializations of the functors of the same name
// in the std lib from C++14 onwards:
//
    struct negate {
        template<class T>
        T operator()(T const &x) const {
            return -x;
        }
    };

    struct plus {
        template<class T, class U>
        T operator()(T const &x, U const &y) const {
            return x + y;
        }
    };

    struct minus {
        template<class T, class U>
        T operator()(T const &x, U const &y) const {
            return x - y;
        }
    };

} // namespace detail

/**
* Returns the zero element for multiplication of polynomials.
*/
template<class T>
Polynomial<T> zero_element(std::multiplies<Polynomial<T> >) {
    return Polynomial<T>();
}

template<class T>
Polynomial<T> identity_element(std::multiplies<Polynomial<T> >) {
    return Polynomial<T>(T(1));
}

/* Calculates a / b and a % b, returning the pair (quotient, remainder) together
* because the same amount of computation yields both.
* This function is not defined for division by zero: user beware.
*/
template<typename T>
std::pair<Polynomial<T>, Polynomial<T> >
quotient_remainder(const Polynomial<T> &dividend, const Polynomial<T> &divisor) {
    assert(divisor);
    if (dividend.size() < divisor.size())
        return std::make_pair(Polynomial<T>(), dividend);
    return detail::division(dividend, divisor);
}


template<class T>
class Polynomial {
public:
    typedef typename std::vector<T>::value_type value_type;
    typedef typename std::vector<T>::size_type size_type;

    Polynomial() {}

    template<class U>
    Polynomial(const U *data, unsigned order)
            : m_data(data, data + order + 1) {
        normalize();
    }

    template<class I>
    Polynomial(I first, I last)
            : m_data(first, last) {
        normalize();
    }

    Polynomial(std::vector<T> &&p) : m_data(std::move(p)) {
        normalize();
    }

    template<class U>
    explicit Polynomial(const U &point, U* = 0) {
        if (point != U(0))
            m_data.push_back(point);
    }

    Polynomial(Polynomial &&p) noexcept
            : m_data(std::move(p.m_data)) {}

    Polynomial(const Polynomial &p)
            : m_data(p.m_data) {}

    template<class U>
    Polynomial(const Polynomial<U> &p) {
        m_data.resize(p.size());
        for (unsigned i = 0; i < p.size(); ++i) {
            m_data[i] = static_cast<T>(p[i]);
        }
    }

    Polynomial(std::initializer_list<T> l) : Polynomial(std::begin(l), std::end(l)) {
    }

    Polynomial &
    operator=(std::initializer_list<T> l) {
        m_data.assign(std::begin(l), std::end(l));
        normalize();
        return *this;
    }

    size_type size() const { return m_data.size(); }

    size_type degree() const {
        if (size() == 0)
            throw std::logic_error("degree() is undefined for the zero polynomial.");
        return m_data.size() - 1;
    }

    value_type &operator[](size_type i) {
        return m_data[i];
    }

    const value_type &operator[](size_type i) const {
        return m_data[i];
    }

    std::vector<T> const &data() const {
        return m_data;
    }

    std::vector<T> &data() {
        return m_data;
    }

    Polynomial &operator=(Polynomial &&p) noexcept {
        m_data = std::move(p.m_data);
        return *this;
    }

    Polynomial &operator=(const Polynomial &p) {
        m_data = p.m_data;
        return *this;
    }

    template<class U>
    Polynomial & operator+=(const U &value) {
        addition(value);
        normalize();
        return *this;
    }

    template<class U>
    Polynomial & operator-=(const U &value) {
        subtraction(value);
        normalize();
        return *this;
    }

    template<class U>
    Polynomial & operator*=(const U &value) {
        multiplication(value);
        normalize();
        return *this;
    }

    template<class U>
    Polynomial & operator/=(const U &value) {
        division(value);
        normalize();
        return *this;
    }

    template<class U>
    Polynomial &
    operator%=(const U & /*value*/) {
        // We can always divide by a scalar, so there is no remainder:
        this->set_zero();
        return *this;
    }

    template<class U>
    Polynomial &operator+=(const Polynomial<U> &value) {
        addition(value);
        normalize();
        return *this;
    }

    template<class U>
    Polynomial &operator-=(const Polynomial<U> &value) {
        subtraction(value);
        normalize();
        return *this;
    }

    template<typename U, typename V>
    void multiply(const Polynomial<U> &a, const Polynomial<V> &b) {
        if (!a || !b) {
            this->set_zero();
            return;
        }
        std::vector<T> prod(a.size() + b.size() - 1, T(0));
        for (unsigned i = 0; i < a.size(); ++i)
            for (unsigned j = 0; j < b.size(); ++j)
                prod[i + j] += a.m_data[i] * b.m_data[j];
        m_data.swap(prod);
    }

    template<class U>
    Polynomial &operator*=(const Polynomial<U> &value) {
        this->multiply(*this, value);
        return *this;
    }

    template<typename U>
    Polynomial &operator/=(const Polynomial<U> &value) {
        *this = quotient_remainder(*this, value).first;
        return *this;
    }

    template<typename U>
    Polynomial &operator%=(const Polynomial<U> &value) {
        *this = quotient_remainder(*this, value).second;
        return *this;
    }

    template<typename U>
    Polynomial &operator>>=(U const &n) {
        assert(n <= m_data.size());
        m_data.erase(m_data.begin(), m_data.begin() + n);
        return *this;
    }

    template<typename U>
    Polynomial &operator<<=(U const &n) {
        m_data.insert(m_data.begin(), n, static_cast<T>(0));
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
        m_data.erase(std::find_if(m_data.rbegin(), m_data.rend(), [](const T &x) -> bool { return x != T(0); }).base(),
                     m_data.end());
    }

private:
    template<class U, class R>
    Polynomial &addition(const U &value, R op) {
        if (m_data.size() == 0)
            m_data.resize(1, 0);
        m_data[0] = op(m_data[0], value);
        return *this;
    }

    template<class U>
    Polynomial &addition(const U &value) {
        return addition(value, detail::plus());
    }

    template<class U>
    Polynomial &subtraction(const U &value) {
        return addition(value, detail::minus());
    }

    template<class U, class R>
    Polynomial &addition(const Polynomial<U> &value, R op) {
        if (m_data.size() < value.size())
            m_data.resize(value.size(), 0);
        for (size_type i = 0; i < value.size(); ++i)
            m_data[i] = op(m_data[i], value[i]);
        return *this;
    }

    template<class U>
    Polynomial &addition(const Polynomial<U> &value) {
        return addition(value, detail::plus());
    }

    template<class U>
    Polynomial &subtraction(const Polynomial<U> &value) {
        return addition(value, detail::minus());
    }

    template<class U>
    Polynomial &multiplication(const U &value) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](const T &x) -> T { return x * value; });
        return *this;
    }

    template<class U>
    Polynomial &division(const U &value) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](const T &x) -> T { return x / value; });
        return *this;
    }

    std::vector<T> m_data;
};


template<class T>
Polynomial<T> operator+(const Polynomial<T> &a, const Polynomial<T> &b) {
    Polynomial<T> result(a);
    result += b;
    return result;
}

template<class T>
Polynomial<T> operator+(Polynomial<T> &&a, const Polynomial<T> &b) {
    a += b;
    return a;
}

template<class T>
Polynomial<T> operator+(const Polynomial<T> &a, Polynomial<T> &&b) {
    b += a;
    return b;
}

template<class T>
Polynomial<T> operator+(Polynomial<T> &&a, Polynomial<T> &&b) {
    a += b;
    return a;
}

template<class T>
Polynomial<T> operator-(const Polynomial<T> &a, const Polynomial<T> &b) {
    Polynomial<T> result(a);
    result -= b;
    return result;
}

template<class T>
Polynomial<T> operator-(Polynomial<T> &&a, const Polynomial<T> &b) {
    a -= b;
    return a;
}

template<class T>
Polynomial<T> operator-(const Polynomial<T> &a, Polynomial<T> &&b) {
    b -= a;
    return -b;
}

template<class T>
Polynomial<T> operator-(Polynomial<T> &&a, Polynomial<T> &&b) {
    a -= b;
    return a;
}

template<class T>
Polynomial<T> operator*(const Polynomial<T> &a, const Polynomial<T> &b) {
    Polynomial<T> result;
    result.multiply(a, b);
    return result;
}

template<class T>
Polynomial<T> operator/(const Polynomial<T> &a, const Polynomial<T> &b) {
    return quotient_remainder(a, b).first;
}

template<class T>
Polynomial<T> operator%(const Polynomial<T> &a, const Polynomial<T> &b) {
    return quotient_remainder(a, b).second;
}

template<class T, class U>
Polynomial<T>
operator+(Polynomial<T> a, const U &b) {
    a += b;
    return a;
}

template<class T, class U>
Polynomial<T>
operator-(Polynomial<T> a, const U &b) {
    a -= b;
    return a;
}

template<class T, class U>
Polynomial<T>
operator*(Polynomial<T> a, const U &b) {
    a *= b;
    return a;
}

template<class T, class U>
Polynomial<T>
operator/(Polynomial<T> a, const U &b) {
    a /= b;
    return a;
}

template<class T, class U>
Polynomial<T>
operator%(const Polynomial<T> &, const U &) {
    return Polynomial<T>();
}

template<class U, class T>
Polynomial<T>
operator+(const U &a, Polynomial<T> b) {
    b += a;
    return b;
}

template<class U, class T>
Polynomial<T>
operator-(const U &a, Polynomial<T> b) {
    b -= a;
    return -b;
}

template<class U, class T>
Polynomial<T>
operator*(const U &a, Polynomial<T> b) {
    b *= a;
    return b;
}

template<class T>
bool operator==(const Polynomial<T> &a, const Polynomial<T> &b) {
    return a.data() == b.data();
}

template<class T>
bool operator!=(const Polynomial<T> &a, const Polynomial<T> &b) {
    return a.data() != b.data();
}

template<typename T, typename U>
Polynomial<T> operator>>(Polynomial<T> a, const U &b) {
    a >>= b;
    return a;
}

template<typename T, typename U>
Polynomial<T> operator<<(Polynomial<T> a, const U &b) {
    a <<= b;
    return a;
}

template<class T>
Polynomial<T> operator-(Polynomial<T> a) {
    std::transform(a.data().begin(), a.data().end(), a.data().begin(), detail::negate());
    return a;
}

template<class charT, class traits, class T>
std::basic_ostream<charT, traits> &operator<<(std::basic_ostream<charT, traits> &os, const Polynomial<T> &poly) {
    os << "{ ";
    for (unsigned i = 0; i < poly.size(); ++i) {
        if (i) os << ", ";
        os << poly[i];
    }
    os << " }";
    return os;
}

#endif //POLYNOMIAL_HPP
