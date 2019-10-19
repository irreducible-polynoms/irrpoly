#ifndef GF_HPP
#define GF_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>

template<uint32_t P = 2>
class gf {
    int_fast64_t v;

    void fix() noexcept;

public:
    static
    gf<P> random() noexcept;

    constexpr
    gf() noexcept;

    gf(int_fast64_t) noexcept;

    int_fast64_t data() const noexcept;

    gf<P> operator+() const noexcept;

    gf<P> operator+(const gf<P> &) const noexcept;

    gf<P> &operator+=(const gf<P> &) noexcept;

    gf<P> &operator++() noexcept;

    gf<P> operator++(int) & noexcept;

    gf<P> operator-() const noexcept;

    gf<P> operator-(const gf<P> &) const noexcept;

    gf<P> &operator-=(const gf<P> &) noexcept;

    gf<P> &operator--() noexcept;

    gf<P> operator--(int) & noexcept;

    gf<P> operator*(const gf<P> &) const noexcept;

    gf<P> &operator*=(const gf<P> &) noexcept;

    [[nodiscard]]
    gf<P> mul_inv() const noexcept(false);

    gf<P> operator/(const gf<P> &) const noexcept(false);

    gf<P> &operator/=(const gf<P> &) noexcept(false);

    [[nodiscard]]
    bool is_zero() const noexcept;

    bool operator==(const gf<P> &) const noexcept;

    bool operator!=(const gf<P> &) const noexcept;

    bool operator>(const gf<P> &) const noexcept;

    bool operator>=(const gf<P> &) const noexcept;

    bool operator<(const gf<P> &) const noexcept;

    bool operator<=(const gf<P> &) const noexcept;

    template<uint32_t Q>
    friend
    std::ostream &operator<<(std::ostream &, const gf<Q> &);

    template<uint32_t Q>
    friend
    std::istream &operator>>(std::istream &, gf<Q> &);
};

template<uint32_t P>
void gf<P>::fix() noexcept {
    v = v % P;
    v = v >= 0 ? v : P + v;
}

template<uint32_t P>
gf<P> gf<P>::random() noexcept {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint_fast64_t> dis(0, P - 1);
    return gf<P>(dis(gen));
}

template<uint32_t P>
constexpr
gf<P>::gf() noexcept : v(0) {}

template<uint32_t P>
gf<P>::gf(const int_fast64_t val) noexcept : v(val) {
    if (P < 2) {
        throw std::domain_error("P must be > 1");
    }
    fix();
}

template<uint32_t P>
int_fast64_t gf<P>::data() const noexcept {
    return v;
}

template<uint32_t P>
gf<P> gf<P>::operator+() const noexcept {
    return gf<P>(+v);
}

template<uint32_t P>
gf<P> &gf<P>::operator+=(const gf<P> &val) noexcept {
    v += val.v;
    fix();
    return *this;
}

template<uint32_t P>
gf<P> gf<P>::operator+(const gf<P> &val) const noexcept {
    return gf<P>(*this) += val;
}

template<uint32_t P>
gf<P> &gf<P>::operator++() noexcept {
    return *this += gf<P>(1);
}

template<uint32_t P>
gf<P> gf<P>::operator++(int) & noexcept {
    gf<P> tmp(*this);
    operator++();
    return tmp;
}

template<uint32_t P>
gf<P> gf<P>::operator-() const noexcept {
    return gf<P>(-v);
}

template<uint32_t P>
gf<P> &gf<P>::operator-=(const gf<P> &val) noexcept {
    v -= val.v;
    fix();
    return *this;
}

template<uint32_t P>
gf<P> gf<P>::operator-(const gf<P> &val) const noexcept {
    return gf<P>(*this) -= val;
}

template<uint32_t P>
gf<P> &gf<P>::operator--() noexcept {
    return *this -= gf<P>(1);
}

template<uint32_t P>
gf<P> gf<P>::operator--(int) & noexcept {
    gf<P> tmp(*this);
    operator++();
    return tmp;
}

template<uint32_t P>
gf<P> &gf<P>::operator*=(const gf<P> &val) noexcept {
    v *= val.v;
    fix();
    return *this;
}

template<uint32_t P>
gf<P> gf<P>::operator*(const gf<P> &val) const noexcept {
    return gf<P>(*this) *= val;
}

template<uint32_t P>
[[nodiscard]]
gf<P> gf<P>::mul_inv() const noexcept(false) {
    static std::array<gf<P>, P> arr{};
    int_fast64_t u0 = P, u1 = 1, u2 = 0, v0 = v, v1 = 0, v2 = 1, w0, w1, w2, q;
    switch (v) {
        case 1:
            return gf<P>(*this);
        case 0:
            throw std::logic_error("multiplicative inverse not exists");
        default:
            switch (arr[v].v) {
                case 1: // marked as irreversible
                    throw std::logic_error("multiplicative inverse not exists");
                case 0: // not calculated
                    while (v0 > 0) {
                        q = u0 / v0;
                        w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
                        u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
                    }
                    if (u0 > 1) {
                        arr[v] = 1; // mark as irreversible
                        throw std::logic_error("multiplicative inverse not exists");
                    }
                    arr[v] = u2; // calculated is inverse for this
                    arr[arr[v].v] = v; // this is inverse for calculated
                default:
                    return arr[v];
            }
    }
}

template<uint32_t P>
gf<P> &gf<P>::operator/=(const gf<P> &val) noexcept(false) {
    if (val.v == 0) { throw std::invalid_argument("division by zero"); }
    return *this *= val.mul_inv();
}

template<uint32_t P>
gf<P> gf<P>::operator/(const gf<P> &val) const noexcept(false) {
    return gf<P>(*this) /= val;
}

template<uint32_t P>
[[nodiscard]]
bool gf<P>::is_zero() const noexcept {
    return v == 0;
}

template<uint32_t P>
bool gf<P>::operator==(const gf<P> &val) const noexcept {
    return v == val.v;
}

template<uint32_t P>
bool gf<P>::operator!=(const gf<P> &val) const noexcept {
    return v != val.v;
}

template<uint32_t P>
bool gf<P>::operator>(const gf<P> &val) const noexcept {
    return v > val.v;
}

template<uint32_t P>
bool gf<P>::operator>=(const gf<P> &val) const noexcept {
    return v >= val.v;
}

template<uint32_t P>
bool gf<P>::operator<(const gf<P> &val) const noexcept {
    return v < val.v;
}

template<uint32_t P>
bool gf<P>::operator<=(const gf<P> &val) const noexcept {
    return v <= val.v;
}

template<uint32_t P>
std::ostream &operator<<(std::ostream &os, const gf<P> &val) {
    return os << val.v;
}

template<uint32_t P>
std::istream &operator>>(std::istream &is, gf<P> &val) {
    return is >> val.v;
}

#endif //GF_HPP
