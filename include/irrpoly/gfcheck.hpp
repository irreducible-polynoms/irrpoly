/**
 * @file    polynomialgf.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#ifndef POLYNOMIALGF_HPP
#define POLYNOMIALGF_HPP

#include <algorithm>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "gfpoly.hpp"
#include "checker.hpp"

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
        gfpoly u0 = m, u1 = gfpoly(m.field(), 1), u2 = gfpoly(m.field(), 0),
               v0 = n, v1 = gfpoly(m.field(), 0), v2 = gfpoly(m.field(), 1),
               w0 = gfpoly(m.field()), w1 = gfpoly(m.field()), w2 = gfpoly(m.field()), q = gfpoly(m.field());
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

        /// Вычисляет производную данного многочлена.
        gfpoly derivative(const gfpoly &val) {
            gfpoly res = val;
            auto i = val.degree();
            for (res[i] = gfn(val.field(), 0); i > 0; --i) {
                res[i - 1] = gfn(val.field(), i) * val[i];
            }
            res.normalize();
            return res;
        }

        /**
         * Вычисляет значение (x^pow - sub) % mod.
         * @param pow степень, в которую требуется возвести x
         * @param mod многочлен, остаток деления на который необходимо найти
         * @param sub вычитаемое, в случае, когда степень многочлена sub меньше
         * степени многочлена mod, можно заменить (x^pow - sub) % mod на
         * (x^pow % mod) - sub без изменения результата, таким образом использование
         * данной функции возможно только в подобной ситуации; это и происходит,
         * поскольку в методе Берлекампа она вызывается всегда с sub = 0,
         * в проверке на примитивность всегда с sub равным константе, при этом
         * mod - как минимум первой степени, поэтому условие выполнено,
         * в методе Рабина sub = x, но при этом mod как минимум второй степени,
         * т.к. все многочлены первой степени неприводимы, что обеспечивает
         * возврат не доходя до вызова данной функции
         */
        [[nodiscard]]
        gfpoly x_pow_mod_sub(uintmax_t pow, const gfpoly &mod, const gfpoly &sub) {
            // Эта функция возврадает по сути следующий редультат:
            // return ((gfpoly({1}) << pow) - sub) % mod;
            // её существование необходимо, т.к. в методе Рабина
            // и в методе проверки на примитивность требуется вычислить
            // её результат для столь больших pow, что x^pow не влезает
            // в оперативную память, данный метод решает эту проблему,
            // но он всё ещё адски медленный
            assert(mod.field() == sub.field());
            const auto n = mod.degree();
            gfpoly xn = gfpoly(mod.field(), 1) << n; // x^n
            gfpoly res(mod.field(), 1);

            uintmax_t d = 0;
            typename gfpoly::size_type tmp;
            for (auto m = res.degree(); pow + m >= n; m = res.degree()) {
                tmp = n - m;
                pow -= tmp;
                res <<= tmp;
                if (res == xn) {
                    // деление в столбик начинается с деления x^n
                    if (d == 0) { d = pow; }
                        // при больших pow деление успевает зациклиться и мы снова придём к x^n
                    else {
                        // находим потерю в степени, произошедшую до этого момента
                        d -= pow;
                        // можем пропустить кучу бесполезных шагов, повторяющих уже сделанные
                        pow %= d;
                        // после этой операции в этот if мы больше не попадём,
                        // т.к. pow < d и деление в столбик не дойдёт до x^n снова,
                        // поэтому d можно не обнулять
                    }
                }
                res %= mod;
            }
            return (res << pow) - sub;
        }

    } // namespace detail

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
        if (n == 0 || (val[0].is_zero() && n > 1)) { return false; }
        if (n == 1) { return true; }

        // функция для построения матрицы берлекампа и вычисления её ранга
        auto berlekampMatrixRank = [](const gfpoly &val) {
            gfpoly tmp(val.field());
            typename gfpoly::size_type i, j, k, l;
            const gfn zer(val.field());
            const auto n = val.degree();
            std::vector<std::vector<gfn>> m(n, std::vector<gfn>(n, zer)); // M = 0
            for (i = 0; i < n; ++i) {
                // M[i,*] = x ^ ip (mod val)
                tmp = detail::x_pow_mod_sub(i * val.field()->base(), val, gfpoly(val.field(), 0));
                for (j = 0, k = tmp.degree(); j <= k; ++j) {
                    m[i][j] += tmp[j];
                }
                m[i][i] -= 1; // M - E
            }

            // приведение матрицы к ступенчатому виду
            bool f;
            gfn mul(val.field());
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
        if (val.is_zero()) { return false; }
        const auto n = val.degree();

        // проверка вырожденных случаев
        if (n == 0 || (val[0].is_zero() && n > 1)) { return false; }
        if (n == 1) { return true; }

        // функция разложения числа на множители
        auto get_list = [](uintmax_t n) {
            std::vector<uintmax_t> list;
            const auto begin = n;
            for (uintmax_t d = 2; d * d <= n; ++d) {
                if (n % d) { continue; }
                list.emplace_back(begin / d);
                while (n % d == 0) { n /= d; }
            }
            if (n != 1) { list.emplace_back(begin / n); }
            return list;
        };

        // шаги 1-2
        auto P = val.field()->base();
        auto list = get_list(n);
        gfpoly tmp(val.field()), x = gfpoly(val.field(), {0, 1});
        for (auto i: list) {
            tmp = detail::x_pow_mod_sub(detail::integer_power(P, i), val, x);
            if (tmp.is_zero() || gcd(val, tmp).degree() > 0) { return false; }
        }

        // шаг 3
        tmp = detail::x_pow_mod_sub(detail::integer_power(P, n), val, x);
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
        if (val.is_zero()) { return false; }
        const auto n = val.degree();

        // проверка вырожденных случаев
        if (n == 0 || (val[0].is_zero() && n > 1)) { return false; }
        if (n == 1) { return true; }

        auto P = val.field()->base();
        gfpoly tmp(val.field()), x = gfpoly(val.field(), {0, 1});
        for (auto m = n / 2, i = 1ul; i <= m; ++i) {
            tmp = detail::x_pow_mod_sub(detail::integer_power(P, i), val, x);
            if (tmp.is_zero() || gcd(val, tmp).degree() > 0) { return false; }
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
        if (val.is_zero()) { return false; }
        const auto n = val.degree();

        // проверка вырожденных случаев
        if (n == 0 || (val[0].is_zero() && n > 1)) { return false; }
        if (n == 1 && val[0] == 0) { return true; } // val = k * x + 0

        // выполняется нормировка, т.к. данный алгоритм справедлив только
        // для многочленов со старшим коэффициентом, равным единице
        // умножение многочлена на число не меняет его примитивность
        const auto poly = val / val[n];

        // ещё один вырожденный случай, на работу с которым алгоритм не рассчитан
        auto P = poly.field()->base();
        if (P == 2 && poly == gfpoly(poly.field(), {1, 1})) { return false; }

        gfn mp = (n % 2) ? -poly[0] : poly[0];

        // функция для разложения (факторизации) числа на множители
        // единица и само число (в случае его простоты) в разложение не входят
        auto factorize = [](uintmax_t n) {
            std::vector<uintmax_t> list;
            const auto begin = n;
            for (uintmax_t d = 2; d * d <= n; ++d) {
                if (n % d) { continue; }
                list.emplace_back(d);
                while (n % d == 0) { n /= d; }
            }
            if (n != 1 && n != begin) { list.emplace_back(n); }
            return list;
        };

        // проверяется выполнение первого условия
        if (P > 2) {
            const auto p = P - 1;
            auto list = (p == 2) ? std::vector<uintmax_t>{2} : factorize(p);
            auto m = list.size() - 1;
            auto tmp = mp;
            for (uint32_t i = 1; i <= p; ++i, tmp *= mp) {
                if (i != p / list[m]) { continue; }
                if (tmp == 1) { return false; }
                if (m == 0) { break; } else { m -= 1; }
            }
        }

        // проверяется выполнение второго условия
        uintmax_t r = (detail::integer_power(P, n) - 1) / (P - 1);
        auto tmp = detail::x_pow_mod_sub(r, val, gfpoly(mp));
        if (!tmp.is_zero()) { return false; }

        // проверяется выполнение третьего условия
        auto list3 = factorize(r);
        const auto m = list3.size();
        for (size_t i = 0; i < m; ++i) {
            tmp = detail::x_pow_mod_sub(r / list3[i], poly, gfpoly(poly.field(), 0));
            if (tmp.is_zero() || tmp.degree() == 0) { return false; }
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
                    case irreducible_method::berlekamp:
                        res.value().irreducible = is_irreducible_berlekamp(poly);
                        break;
                    case irreducible_method::rabin:
                        res.value().irreducible = is_irreducible_rabin(poly);
                        break;
                    case irreducible_method::benor:
                        res.value().irreducible = is_irreducible_benor(poly);
                        break;
                    default: ; // irreducible_method::nil
                }
                switch (prim_meth) {
                    case primitive_method::definition:
                        res.value().primitive = res.value().irreducible ? is_primitive_definition(poly) : false;
                        break;
                    default: ; // primitive_method::nil
                }
            };
        }

    } // namespace multithread

} // namespace irrpoly

#endif //POLYNOMIALGF_HPP
