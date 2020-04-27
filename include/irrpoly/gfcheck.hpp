/**
 * @file    polynomialgf.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include "checker.hpp"
#include "gfpoly.hpp"

#include <algorithm>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace irrpoly {

/**
 * Расширенный алгоритм Евклида для поиска наибольшего общего делителя
 * (greatest common divisor) двух многочленов. Реализация сделана на основе
 * кода из библиотеки Boost 1.71.0.
 */
gfpoly gcd(gfpoly m, gfpoly n) {
    assert(m.field() == n.field());
    if (m.is_zero() || n.is_zero()) {
        throw std::domain_error("arguments must be strictly positive");
    }
    if (m.degree() < n.degree()) {
        std::swap(m, n);
    }
    gfpoly u0 = m, u1 = gfpoly(m.field(), 1), u2 = gfpoly(m.field()),
        v0 = n, v1 = gfpoly(m.field()), v2 = gfpoly(m.field(), 1),
        w0 = gfpoly(m.field()), w1 = gfpoly(m.field()), w2 = gfpoly(m.field()),
        q = gfpoly(m.field());
    while (!v0.is_zero()) {
        q = u0 / v0;
        w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
        u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
    }
    return u0;
}

namespace detail {

/// Быстрое возведение в степень, взято из библиотеки Boost.
template<class T, class N>
T integer_power(T t, N n) {
    switch (n) {
    case 0:return static_cast<T>(1u);
    case 1:return t;
    case 2:return t * t;
    case 3:return t * t * t;
    default:T result = integer_power(t, n / 2);
        return (n & 1u) ? result * result * t : result * result;
    }
}

/// Вычисляет производную данного многочлена.
gfpoly derivative(const gfpoly &val) {
    assert(val && val.degree());
    std::vector<uintmax_t> res(val.size() - 1, 0);
    for (uintmax_t i = 1; i < val.size(); ++i) {
        res[i - 1] = (i * val[i]).value();
    }
    return gfpoly(val.field(), res);
}

/**
 * Вычисляет значение (x^pow) % mod.
 * Позволяет экономить память за счёт потери в скорости.
 * Принцип работы: деление в столбик начинается с деления x^n,
 * при больших pow деление успевает зациклиться и мы снова придём к x^n.
 * @param pow степень, в которую требуется возвести x
 * @param mod многочлен, остаток деления на который необходимо найти
 */
[[nodiscard]]
gfpoly x_pow_mod(uintmax_t pow, const gfpoly &mod) {
    const auto n = mod.degree();
    gfpoly xn = gfpoly(mod.field(), 1) << n; // x^n
    gfpoly res(mod.field(), 1);

    uintmax_t d = 0;
    uintmax_t tmp;
    for (auto m = res.degree(); pow + m >= n; m = res.degree()) {
        tmp = n - m, pow -= tmp, res <<= tmp;
        d = (res == xn) ? (d) ? (d -= pow, pow %= d, 0) : pow : d;
        res %= mod;
    }
    return res << pow;
}

} // namespace

/**
 * Алгоритм Берлекампа проверки многочлена на неприводимость в поле GF[P].
 * Первый шаг - вычисление производной данного многочлена. Если производная
 * равна нулю, то многочлен является степенью какого-то другого многочлена,
 * то есть он приводим.
 * Второй шаг - поиск общих множителей многочлена и его производной.
 * Если общие множители (многочлены, а не числа) есть, т.е. многочлены
 * не взаимно просты, то val делится на них, т.е. он не неприводим.
 * Третий шаг - простоение матрицы Берлекампа и вычисление её ранга.
 * Строится матрица M[nxn], где строки - коэффициенты многочлена x^(iP) (mod val),
 * где P - основание поля галуа, 0 < i < n, val - текущий многочлен над полем GF[P].
 * Подробное описание и пример расчёта можно найти в статье
 * "A Formalization of Berlekamp’s Factorization Algorithm" по ссылке
 * http://www21.in.tum.de/~nipkow/Isabelle2016/Isabelle2016_6.pdf (стр. 3-4).
 * Из матрицы M вычитается единичная матрица и получается матрица Берлекампа.
 * Если ранг матрицы Берлекампа равен степени многочлена минус 1,
 * то многочлен неприводим. Для вычисления ранга используется приведение
 * матрицы к ступенчатому виду и подсчёт числа ступеней в ней.
 * Кроме того, все многочлены первой степени неприводимы в любом поле.
 */
bool is_irreducible_berlekamp(const gfpoly &val) {
    if (val.is_zero()) {
        return false;
    }
    const auto n = val.degree();

    // проверка вырожденных случаев
    if (n == 0 || (val[0].is_zero() && n > 1)) {
        return false;
    }
    if (n == 1) {
        return true;
    }

    // функция для построения матрицы берлекампа и вычисления её ранга
    auto berlekampMatrixRank = [](const gfpoly &val) {
        gfpoly poly(val.field());
        uintmax_t i, j, k, l;
        const auto n = val.degree();
        std::vector<std::vector<gfn>> m(n, std::vector<gfn>(n, gfn(val.field()))); // M = 0
        for (i = 0; i < n; ++i) {
            // M[i,*] = x ^ ip (mod val)
            poly = detail::x_pow_mod(i * val.base(), val);
            for (j = 0, k = poly.degree(); j <= k; ++j) {
                m[i][j] += poly[j];
            }
            m[i][i] -= 1; // M - E
        }

        // приведение матрицы к ступенчатому виду
        bool f;
        gfn num(val.field());
        for (i = k = 0; i < n && k < n; ++k) {
            f = !m[i][k].is_zero();
            for (j = i + 1; j < n; ++j) {
                if (!m[j][k].is_zero()) {
                    if (f) {
                        num = m[j][k] / m[i][k];
                        m[j][k].set_zero();
                        for (l = k + 1; l < n; ++l) {
                            m[j][l] -= m[i][l] * num;
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

    // алгоритм Берлекампа
    auto d = detail::derivative(val);
    return !d.is_zero() && gcd(val, d).degree() == 0 &&
        berlekampMatrixRank(val) == val.degree() - 1;
}

/**
 * Алгоритм Рабина проверки многочлена на неприводимость в поле GF[P].
 * Для использования алгоритма предварительно строится список простых множителей
 * для числа n - степени многочлена. Вместо множителей в список добавляется отношение
 * n_i = n / d_i, где d_i - простой делитель (divisor) n. Если n - простое, то
 * список состоит из одной единицы. Далее выполняется несколько шагов:
 * 1. строится многочлен temp =  x^(P ^ n_i) - x (mod val)
 * 2. находится наиболиший общий дилитель temp и val, если НОД не равен константе,
 * отличной от нуля, то это многочлен, и val, очевидно, делится на него, а значит приводим.
 * 3. если условия 1-2 выполнены для всех n_i, проверяем их для n.
 * Если для n получаем результат 0, то val неприводим.
 * Кроме того, все многочлены первой степени неприводимы в любом поле.
 * Подробную информацию по алгоритму можно найти здесь:
 * https://www.fing.edu.uy/inco/pedeciba/bibliote/reptec/TR0116.pdf
 * или на странице Википедии в разделе Rabin's test of irreducibility
 * https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
 */
bool is_irreducible_rabin(const gfpoly &val) {
    if (val.is_zero()) {
        return false;
    }
    const auto n = val.degree();

    // проверка вырожденных случаев
    if (n == 0 || (val[0].is_zero() && n > 1)) {
        return false;
    }
    if (n == 1) {
        return true;
    }

    // функция разложения числа на множители
    auto get_list = [](uintmax_t n) {
        std::vector<uintmax_t> list;
        const auto begin = n;
        for (uintmax_t d = 2; d * d <= n; ++d) {
            if (n % d) {
                continue;
            }
            list.emplace_back(begin / d);
            while (n % d == 0) {
                n /= d;
            }
        }
        if (n != 1) {
            list.emplace_back(begin / n);
        }
        return list;
    };

    // шаги 1-2
    auto P = val.base();
    auto list = get_list(n);
    gfpoly tmp(val.field()), x = gfpoly(val.field(), {0, 1});
    for (auto i: list) {
        tmp = detail::x_pow_mod(detail::integer_power(P, i), val) - x;
        if (tmp.is_zero() || gcd(val, tmp).degree() > 0) {
            return false;
        }
    }

    // шаг 3
    tmp = detail::x_pow_mod(detail::integer_power(P, n), val) - x;
    return tmp.is_zero();
}

/**
 * Алгоритм Бен-Ора проверки многочлена на неприводимость в поле GF[P].
 * Для всех i от 1 до m/2 включительно, где m - степень проверяемого многочлена:
 * 1. строится многочлен temp =  x^(P ^ i) - x (mod val)
 * 2. находится наиболиший общий дилитель temp и val, если НОД не равен константе,
 * отличной от нуля, то это многочлен, и val, очевидно, делится на него, а значит приводим.
 * После завершения шагов 1-2 для всех i получаем, что многочлен неприводим.
 * Кроме того, все многочлены первой степени неприводимы в любом поле.
 */
bool is_irreducible_benor(const gfpoly &val) {
    if (val.is_zero()) {
        return false;
    }
    const auto n = val.degree();

    // проверка вырожденных случаев
    if (n == 0 || (val[0].is_zero() && n > 1)) {
        return false;
    }
    if (n == 1) {
        return true;
    }

    auto P = val.base();
    gfpoly tmp(val.field()), x = gfpoly(val.field(), {0, 1});
    for (uintmax_t m = n / 2, i = 1; i <= m; ++i) {
        tmp = detail::x_pow_mod(detail::integer_power(P, i), val) - x;
        if (tmp.is_zero() || gcd(val, tmp).degree() > 0) {
            return false;
        }
    }
    return true;
}

/**
 * Алгоритм проверки многочлена на примитивность по определению. Многочлен является
 * примитивным над полем GF[P], если выполнены три условия:
 * 1. элемент mp = (-1)^n * val[0] является примитивным элементом поля GF[P^n], т.е.
 * k^((p-1) / q) != 1 для каждого q - простого множителя P-1
 * данный пункт не применим для P = 2 по объективным причинам
 * 2. x^r = k (mod val), где r = (p^n - 1) / (p - 1)
 * 3. deg[x^(r / q) (mod val)] > 0 для каждого 1 < q < r - простого множителя r
 * Кроме того, многочлен x является примитивным для любого поля GF[P].
 * Подробную информацию по алгоритму можно найти здесь
 * https://www.ams.org/journals/mcom/1992-59-200/S0025-5718-1992-1134730-7/S0025-5718-1992-1134730-7.pdf
 * Возможные пути параллелизации данного алгоритма приведены в статье
 * https://www.researchgate.net/publication/329358609_Parallelization_of_Algorithm_for_Primitive_Polynomials_Generation_in_Extended_Galois_Field_pm
 */
bool is_primitive_definition(const gfpoly &val) {
    if (val.is_zero()) {
        return false;
    }
    const auto n = val.degree();

    // проверка вырожденных случаев
    if (n == 0 || (val[0].is_zero() && n > 1)) {
        return false;
    }
    if (n == 1 && val[0] == 0) {
        return true;
    } // val = k * x + 0

    // выполняется нормировка, т.к. данный алгоритм справедлив только
    // для многочленов со старшим коэффициентом, равным единице
    // умножение многочлена на число не меняет его примитивность
    const auto poly = val / val[n];

    // ещё один вырожденный случай, на работу с которым алгоритм не рассчитан
    auto P = poly.base();
    if (P == 2 && poly == gfpoly(poly.field(), {1, 1})) {
        return false;
    }

    gfn mp = (n % 2) ? -poly[0] : poly[0];

    // функция для разложения (факторизации) числа на множители
    // единица и само число (в случае его простоты) в разложение не входят
    auto factorize = [](uintmax_t n) {
        std::vector<uintmax_t> list;
        const auto begin = n;
        for (uintmax_t d = 2; d * d <= n; ++d) {
            if (n % d) {
                continue;
            }
            list.emplace_back(d);
            while (n % d == 0) {
                n /= d;
            }
        }
        if (n != 1 && n != begin) {
            list.emplace_back(n);
        }
        return list;
    };

    // проверяется выполнение первого условия
    if (P > 2) {
        const auto p = P - 1;
        auto list = (p == 2) ? std::vector<uintmax_t>{2} : factorize(p);
        auto m = list.size() - 1;
        auto tmp = mp;
        for (uint32_t i = 1; i <= p; ++i, tmp *= mp) {
            if (i != p / list[m]) {
                continue;
            }
            if (tmp == 1) {
                return false;
            }
            if (m == 0) {
                break;
            } else {
                m -= 1;
            }
        }
    }

    // проверяется выполнение второго условия
    uintmax_t r = (detail::integer_power(P, n) - 1) / (P - 1);
    auto tmp = detail::x_pow_mod(r, val) - mp;
    if (!tmp.is_zero()) {
        return false;
    }

    // проверяется выполнение третьего условия
    auto list3 = factorize(r);
    const auto m = list3.size();
    for (size_t i = 0; i < m; ++i) {
        tmp = detail::x_pow_mod(r / list3[i], poly);
        if (tmp.is_zero() || tmp.degree() == 0) {
            return false;
        }
    }

    // если все условия выполнены - многочлен примитивен
    return true;
}

namespace multithread {

/// Структура, представляющая результаты проверки многочлена.
struct result_type {
    bool irreducible;
    bool primitive;
};

/// Доступные методы проверки нанеприводимость.
enum class irreducible_method {
    nil, ///< не проверять
    berlekamp, ///< алгоритм Берлекампа
    rabin, ///< алгоритм Рабина
    benor, ///< алгоритм Бен-Ора
};

/// Доступные методы проверки примитивность.
enum class primitive_method {
    nil, ///< не проверять
    definition, ///< проверка по определению
};

using polychecker = checker<gfpoly, result_type>;

/// Формируется универсальная функция проверки многочленов.
typename checker<gfpoly, result_type>::check_func make_check_func(
    irreducible_method irr_meth, primitive_method prim_meth) {
    return [=](const gfpoly &poly, std::optional<result_type> &res) {
        // в случае, когда проверка не выполняется устанавливается результат true
        res.emplace(result_type{true, true});
        switch (irr_meth) {
        case irreducible_method::berlekamp:res.value().irreducible = is_irreducible_berlekamp(poly);
            break;
        case irreducible_method::rabin:res.value().irreducible = is_irreducible_rabin(poly);
            break;
        case irreducible_method::benor:res.value().irreducible = is_irreducible_benor(poly);
            break;
        default:; // irreducible_method::nil
        }
        switch (prim_meth) {
        case primitive_method::definition:
            res.value().primitive = res.value().irreducible ? is_primitive_definition(poly) : false;
            break;
        default:; // primitive_method::nil
        }
    };
}

} // namespace multithread

} // namespace irrpoly
