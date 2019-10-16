#ifndef GF_HPP
#define GF_HPP

#include <cstdint>
#include <iostream>
#include <stdexcept>

#include "boost/integer/mod_inverse.hpp"

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

    GF<P> MulInv() const noexcept(false);
    
    GF<P> operator/(const GF<P> &) const noexcept(false);

    GF<P> &operator/=(const GF<P> &) noexcept(false);

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
GF<P> GF<P>::MulInv() const noexcept(false) {
    int64_t tmp = boost::integer::mod_inverse<int64_t>(v, P);
    if (tmp == 0) { throw std::logic_error("multiplicative inverse not exists"); }
    return tmp;
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
