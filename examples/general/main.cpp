#include <iostream>
#include <thread>

#include "irrpoly/polynomialgf.hpp"

using namespace irrpoly;

/**
 * Шаблон функции, применяемой для генерации многочленов с требуемыми характеристиками.
 * В функции генерируются случайные многочлены требуемой степени и выполняется их проверка
 * на неприводимость и примитивность. Для осмысленного применения рекоммендуется перебирать
 * многочлены последовательно, а не выбирать каждый раз случайный (во избежание повторов).
 * @param num число многочленов с требуемыми характеристиками, которые необходимо получить
 * @param degree степень многочленов, которые требуется проверять
 * @param irr_meth метод проверки на неприводимость
 * @param prim_meth метод проверки на примитивность
 * @param threads_num число потоков, выполняющих проверку (число физических ядер минус один)
 * @return вектор сгенерированных многочленов с требуемыми характеристиками
 */
template<uint32_t P>
[[nodiscard]]
std::vector<polynomialgf<P>> generate_irreducible(
        const typename std::vector<polynomialgf<P>>::size_type num,
        const typename polynomialgf<P>::size_type degree,
        const typename multithread::irreducible_method irr_meth,
        const typename multithread::primitive_method prim_meth,
        const unsigned threads_num
) {
    // возвращаемое значение
    std::vector<polynomialgf<P>> arr;
    arr.reserve(num);

    // создаём всё необходимое для многопоточности
    multithread::polychecker<P> ch(threads_num);

    auto input = [&]() -> polynomialgf<P> {
        return random < P > (degree);
    };

    auto check = multithread::make_check_func<P>(irr_meth, prim_meth);

    auto callback = [&](const polynomialgf<P> &poly, const typename multithread::result_type &res) -> bool {
        if (res.irreducible) {
            arr.emplace_back(poly);
            return true;
        }
        return false;
    };

    ch.check(num, input, check, callback);

    return arr;
}

int main() {
    // основание поля Галуа GF[P]
    const uint32_t P = 2;
    // число многочленов, которые требуется найти
    const typename std::vector<polynomialgf<P>>::size_type num = 3;
    // степень искомых многочленов
    const typename polynomialgf<P>::size_type degree = 5;
    // какой из методов проверки на неприводимость хотим использовать
    // возможные варианты - отсутствие приверки (nil), Берлекампа (berlekamp) и Рабин (rabin)
    const typename multithread::irreducible_method irr_meth = multithread::irreducible_method::berlekamp;
    // какой из методов проверки на примитивность хотим использовать
    // возможные варианты - отсутствие приверки (nil) и проверка по определению (definition)
    const typename multithread::primitive_method prim_meth = multithread::primitive_method::nil;
    // число потоков, выполняюих проверку многочленов на соответствие заданной характеристике
    // должно быть равно числу физических ядер минус один
    const unsigned threads_num = std::thread::hardware_concurrency() - 1;

    // генерируем многочлены и выводим результат
    auto poly = generate_irreducible<P>(num, degree, irr_meth, prim_meth, threads_num);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }

    return 0;
}
