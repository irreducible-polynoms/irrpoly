#ifndef POLYNOMIALGF_HPP
#define POLYNOMIALGF_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "gf.hpp"
#include "polynomial.hpp"

template<uint32_t P = 2>
using polynomialgf = polynomial<gf<P>>;

template<uint32_t P>
polynomialgf<P> derivative(const polynomialgf<P> &val) {
    polynomialgf<P> res = val;
    auto i = val.degree();
    for (res[i] = 0; i > 0; --i) {
        res[i - 1] = gf<P>(i) * val[i];
    }
    res.normalize();
    return res;
}

template<uint32_t P>
polynomialgf<P> gcd(polynomialgf<P> m, polynomialgf<P> n) {
    if (m.is_zero() || n.is_zero()) {
        throw std::domain_error("arguments must be strictly positive");
    }
    if (m.degree() < n.degree()) {
        std::swap(m, n);
    }
    polynomialgf<P> u0 = m, u1 = polynomialgf<P>({1}), u2 = polynomialgf<P>({0}),
            v0 = n, v1 = polynomialgf<P>({0}), v2 = polynomialgf<P>({1}),
            w0, w1, w2, q;
    while (!v0.is_zero()) {
        q = u0 / v0;
        w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
        u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
    }
    return u0;
}

template<uint32_t P>
bool is_irreducible_berlekamp(const polynomialgf<P> &val) {
    const auto n = val.degree();
    if (val.is_zero() || n == 0 || (val[0].is_zero() && n > 1)) { return false; }
    if (n == 1) { return true; }
    auto berlekampMatrixRank = [](const polynomialgf<P> &val) {
        polynomialgf<P> tmp;
        typename polynomialgf<P>::size_type i, j, k, l;
        const gf<P> zer = 0;
        const auto n = val.degree();
        std::vector<std::vector<gf<P>>> m(n, std::vector<gf<P>>(n, zer)); // berlekamp matrix
        for (i = 0; i < n; ++i) {
            tmp = (polynomialgf<P>({1}) << i * P) % val; // temp = x ^ ip (mod val)
            for (j = 0, k = tmp.degree(); j <= k; ++j) {
                m[i][j] += tmp[j];
            }
            m[i][i] -= 1; // m -= E
        }

        // reduction of a matrix to a step form with calculation of its rank
        bool f;
        gf<P> mul;
        for (i = k = 0; i < n && k < n; ++k) {
            f = !m[i][k].is_zero();
            for (j = i + 1; j < n; ++j) {
                if (!m[j][k].is_zero()) {
                    if (f) {
                        mul = m[i][k].mul_inv() * m[j][k];
                        m[j][k] = zer;
                        for (l = k + 1; l < n; ++l) {
                            m[j][l] -= m[i][l] * mul;
                        }
                    } else {
                        for (l = k; l < n; ++l) {
                            std::swap(m[i][l], m[j][l]);
                        }
                        f = true;
                    }
                }
            }
            i += f;
        }
        return i;
    };
    auto d = derivative(val);
    return !d.is_zero() && gcd(val, d).degree() == 0 &&
           berlekampMatrixRank(val) == val.degree() - 1;
}

template<uint32_t P>
polynomialgf<P> random(typename polynomialgf<P>::size_type degree) {
    std::vector<gf<P>> data(degree + 1);
    for (auto &d : data) { d = gf<P>::random(); }
    while (data[degree].is_zero()) { data[degree] = gf<P>::random(); }
    return data;
}

template<uint32_t P>
bool is_primitive(const polynomialgf<P> &val) {
    const auto n = val.degree();
    if (val.is_zero() || n == 0 || (val[0].is_zero() && n > 1)) { return false; }
    const auto poly = val / val[n];
    if (P == 2 && poly == polynomialgf<P>({ 1,1 })) { return false; }
    auto mp = (n % 2) ? -poly[0] : poly[0];

    auto factorize = [](uint64_t n) {
        std::vector<uint64_t> list;
        const auto begin = n;
        for (uint64_t d = 2; d * d <= n; ++d) {
            if (n % d) { continue; }
            list.emplace_back(d);
            while (n % d == 0) { n /= d; }
        }
        if (n != 1 && n != begin) { list.emplace_back(n); }
        return list;
    };

    // step 1
    if (P > 2) {
        const auto p = P - 1;
        auto list = (p == 2) ? std::vector<uint64_t>{2} : factorize(p);
        auto m = list.size() - 1;
        auto tmp = mp;
        for (uint32_t i = 1; i <= p; ++i, tmp *= mp) {
            if (i != p / list[m]) { continue; }
            if (tmp.data() == 1) { return false; }
            if (m == 0) { break; } else { m -= 1; }
        }
    }

    // step 2
    uint_fast64_t r = 1;
    for (auto i = n; i > 0; --i) { r *= P; }
    r = (r - 1) / (P - 1);
    auto tmp = (polynomialgf<P>({ 1 }) << r) - polynomialgf<P>({ mp });
    if (!(tmp % val).is_zero()) { return false; }

    // step 3
    auto list3 = factorize(r);
    const auto m = list3.size();
    for (size_t i = 0; i < m; ++i) {
        tmp = (polynomialgf<P>({ 1 }) << (r / list3[i])) % poly;
        if (tmp.is_zero() || tmp.degree() == 0) { return false; }
    }

    return true;
}

template<uint32_t P>
bool is_irreducible_rabin(const polynomialgf<P> &val) {
    auto get_list = [](uint64_t n) {
        std::vector<uint64_t> list;
        const auto begin = n;
        for (int64_t d = 2; d * d <= n; ++d) {
            if (n % d) { continue; }
            list.emplace_back(begin / d);
            while (n % d == 0) { n /= d; }
        }
        if (n != 1) { list.emplace_back(begin / n); }
        return list;
    };

    const auto n = val.degree();
    auto list = get_list(n);
    polynomialgf<P> tmp;
    for (auto i: list) {
        tmp = polynomialgf<P>({1}) << pow(P, i);
        tmp -= polynomialgf<P>({1}) << 1;
        tmp %= val;
        if (tmp.is_zero() || gcd(val, tmp).degree() > 0) { return false; }
    }
    tmp = polynomialgf<P>({1}) << pow(P, n);
    tmp -= polynomialgf<P>({1}) << 1;
    tmp %= val;
    return tmp.is_zero();
}

#endif //POLYNOMIALGF_HPP
