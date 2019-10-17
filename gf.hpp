#ifndef GF_HPP
#define GF_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>

template<uint32_t P = 2>
class GF {
    int64_t v;

    void fix() noexcept;

public:
    GF(int64_t = 0) noexcept;

    GF<P> operator+() const noexcept;

    GF<P> operator+(const GF<P> &) const noexcept;

    GF<P> &operator+=(const GF<P> &) noexcept;

    GF<P> &operator++() noexcept;

    GF<P> operator++(int) & noexcept;

    GF<P> operator-() const noexcept;

    GF<P> operator-(const GF<P> &) const noexcept;

    GF<P> &operator-=(const GF<P> &) noexcept;

    GF<P> &operator--() noexcept;

    GF<P> operator--(int) & noexcept;

    GF<P> operator*(const GF<P> &) const noexcept;

    GF<P> &operator*=(const GF<P> &) noexcept;

    [[nodiscard]]
    GF<P> MulInv() const noexcept(false);

    GF<P> operator/(const GF<P> &) const noexcept(false);

    GF<P> &operator/=(const GF<P> &) noexcept(false);

    [[nodiscard]]
    bool IsZero() const noexcept;

    bool operator==(const GF<P> &) const noexcept;

    bool operator!=(const GF<P> &) const noexcept;

    bool operator>(const GF<P> &) const noexcept;

    bool operator>=(const GF<P> &) const noexcept;

    bool operator<(const GF<P> &) const noexcept;

    bool operator<=(const GF<P> &) const noexcept;

    template<uint32_t Q>
    friend
    std::ostream &operator<<(std::ostream &, const GF<Q> &);

    template<uint32_t Q>
    friend
    std::istream &operator>>(std::istream &, GF<Q> &);
};

template<uint32_t P>
void GF<P>::fix() noexcept {
    v = v % P;
    v = v >= 0 ? v : P + v;
}

template<uint32_t P>
GF<P>::GF(const int64_t val) noexcept : v(val) {
    if (P < 2) {
        throw std::domain_error("P must be > 1");
    }
    fix();
}

template<uint32_t P>
GF<P> GF<P>::operator+() const noexcept {
    return GF<P>(+v);
}

template<uint32_t P>
GF<P> &GF<P>::operator+=(const GF<P> &val) noexcept {
    v += val.v;
    fix();
    return *this;
}

template<uint32_t P>
GF<P> GF<P>::operator+(const GF<P> &val) const noexcept {
    return GF<P>(*this) += val;
}

template<uint32_t P>
GF<P> &GF<P>::operator++() noexcept {
    return *this += GF<P>(1);
}

template<uint32_t P>
GF<P> GF<P>::operator++(int) & noexcept {
    GF<P> tmp(*this);
    operator++();
    return tmp;
}

template<uint32_t P>
GF<P> GF<P>::operator-() const noexcept {
    return GF<P>(-v);
}

template<uint32_t P>
GF<P> &GF<P>::operator-=(const GF<P> &val) noexcept {
    v -= val.v;
    fix();
    return *this;
}

template<uint32_t P>
GF<P> GF<P>::operator-(const GF<P> &val) const noexcept {
    return GF<P>(*this) -= val;
}

template<uint32_t P>
GF<P> &GF<P>::operator--() noexcept {
    return *this -= GF<P>(1);
}

template<uint32_t P>
GF<P> GF<P>::operator--(int) & noexcept {
    GF<P> tmp(*this);
    operator++();
    return tmp;
}

template<uint32_t P>
GF<P> &GF<P>::operator*=(const GF<P> &val) noexcept {
    v *= val.v;
    fix();
    return *this;
}

template<uint32_t P>
GF<P> GF<P>::operator*(const GF<P> &val) const noexcept {
    return GF<P>(*this) *= val;
}

template<uint32_t P>
[[nodiscard]]
GF<P> GF<P>::MulInv() const noexcept(false) {
    if (v == 0) { return GF<P>(*this); }
    int64_t u0 = P, u1 = 1, u2 = 0, v0 = v, v1 = 0, v2 = 1, w0, w1, w2, q;
    while (v0 > 0) {
        q = u0 / v0;
        w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
        u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
    }
    if (u0 > 1) {
        throw std::logic_error("multiplicative inverse not exists");
    }
    return GF<P>(u2);
}

template<uint32_t P>
GF<P> &GF<P>::operator/=(const GF<P> &val) noexcept(false) {
    if (val.v == 0) { throw std::invalid_argument("division by zero"); }
    return *this *= val.MulInv();
}

template<uint32_t P>
GF<P> GF<P>::operator/(const GF<P> &val) const noexcept(false) {
    return GF<P>(*this) /= val;
}

template<uint32_t P>
[[nodiscard]]
bool GF<P>::IsZero() const noexcept {
    return v == 0;
}

template<uint32_t P>
bool GF<P>::operator==(const GF<P> &val) const noexcept {
    return v == val.v;
}

template<uint32_t P>
bool GF<P>::operator!=(const GF<P> &val) const noexcept {
    return v != val.v;
}

template<uint32_t P>
bool GF<P>::operator>(const GF<P> &val) const noexcept {
    return v > val.v;
}

template<uint32_t P>
bool GF<P>::operator>=(const GF<P> &val) const noexcept {
    return v >= val.v;
}

template<uint32_t P>
bool GF<P>::operator<(const GF<P> &val) const noexcept {
    return v < val.v;
}

template<uint32_t P>
bool GF<P>::operator<=(const GF<P> &val) const noexcept {
    return v <= val.v;
}

template<uint32_t P>
std::ostream &operator<<(std::ostream &os, const GF<P> &val) {
    return os << val.v;
}

template<uint32_t P>
std::istream &operator>>(std::istream &is, GF<P> &val) {
    return is >> val.v;
}

#endif //GF_HPP
