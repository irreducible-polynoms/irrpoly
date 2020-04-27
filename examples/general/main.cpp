// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include <irrpoly.h>

#include <iostream>
#include <thread>

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
[[nodiscard]]
auto generate_irreducible(
    const uintmax_t base,
    const uintmax_t num,
    const uintmax_t degree,
    const typename multithread::irreducible_method irr_meth,
    const typename multithread::primitive_method prim_meth,
    const unsigned threads_num
) -> std::vector<gfpoly> {
    // возвращаемое значение
    std::vector<gfpoly> arr;
    arr.reserve(num);

    // создаём всё необходимое для многопоточности
    multithread::polychecker ch(threads_num);

    auto field = make_gf(base);
    auto input = [&]() -> gfpoly {
        return gfpoly::random(field, degree);
    };

    auto check = multithread::make_check_func(irr_meth, prim_meth);

    auto callback = [&](const gfpoly &poly, const typename multithread::result_type &res) -> bool {
        if (res.irreducible) {
            arr.emplace_back(poly);
            return true;
        }
        return false;
    };

    ch.check(num, input, check, callback);

    return arr;
}

auto main() -> int {
    // основание поля Галуа
    const uintmax_t base = 2;
    // число многочленов, которые требуется найти
    const uintmax_t num = 3;
    // степень искомых многочленов
    const uintmax_t degree = 5;
    // какой из методов проверки на неприводимость хотим использовать
    // возможные варианты - отсутствие приверки (nil), Берлекампа (berlekamp), Рабина (rabin) и Бенора (benor)
    const typename multithread::irreducible_method irr_meth = multithread::irreducible_method::benor;
    // какой из методов проверки на примитивность хотим использовать
    // возможные варианты - отсутствие приверки (nil) и проверка по определению (definition)
    const typename multithread::primitive_method prim_meth = multithread::primitive_method::nil;
    // число потоков, выполняюих проверку многочленов на соответствие заданной характеристике
    // должно быть равно числу физических ядер минус один
    const unsigned threads_num = std::thread::hardware_concurrency();

    // генерируем многочлены и выводим результат
    auto poly = generate_irreducible(base, num, degree, irr_meth, prim_meth, threads_num);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }

    return 0;
}
