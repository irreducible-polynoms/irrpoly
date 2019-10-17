#ifndef POLYNOMIALGF_HPP
#define POLYNOMIALGF_HPP

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "gf.hpp"
#include "polynomial.hpp"

template<uint32_t P = 2>
using PolynomialGF = Polynomial<GF<P>>;

template<uint32_t P>
PolynomialGF<P> derivative(const PolynomialGF<P> &val) {
    PolynomialGF<P> res = val;
    auto i = val.degree();
    for (res[i] = 0; i > 0; --i) {
        res[i - 1] = GF<P>(i) * val[i];
    }
    res.normalize();
    return res;
}

template<uint32_t P>
PolynomialGF<P> gcd(PolynomialGF<P> m, PolynomialGF<P> n) {
    if (m.is_zero() || n.is_zero()) {
        throw std::domain_error("arguments must be strictly positive");
    }
    if (m.degree() < n.degree()) {
        std::swap(m, n);
    }
    PolynomialGF<P> u0 = m, u1 = PolynomialGF<P>({1}), u2 = PolynomialGF<P>({0}),
            v0 = n, v1 = PolynomialGF<P>({0}), v2 = PolynomialGF<P>({1}),
            w0, w1, w2, q;
    while (!v0.is_zero()) {
        q = u0 / v0;
        w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
        u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
    }
    return u0;
}

template<uint32_t P>
bool isIrreducible(const PolynomialGF<P> &val) {
    if (val.degree() == 0) { return true; }
    if (val[0].IsZero()) { return false; }
    auto berlekampMatrixRank = [](const PolynomialGF<P> &val) {
        PolynomialGF<P> tmp;
        typename PolynomialGF<P>::size_type i, j, k, l;
        const GF<P> zer = 0;
        const auto n = val.degree();
        std::vector<std::vector<GF<P>>> m(n, std::vector<GF<P>>(n, zer)); // berlekamp matrix
        for (i = 0; i < n; ++i) {
            tmp = (PolynomialGF<P>({1}) << i * P) % val; // temp = x ^ ip (mod val)
            for (j = 0, k = tmp.degree(); j <= k; ++j) {
                m[i][j] += tmp[j];
            }
            m[i][i] -= 1; // m -= E
        }

        // reduction of a matrix to a step form with calculation of its rank
        bool f;
        GF<P> mul;
        for (i = k = 0; i < n && k < n; ++k) {
            f = !m[i][k].IsZero();
            for (j = i + 1; j < n; ++j) {
                if (!m[j][k].IsZero()) {
                    if (f) {
                        mul = m[i][k].MulInv() * m[j][k];
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

#endif //POLYNOMIALGF_HPP
