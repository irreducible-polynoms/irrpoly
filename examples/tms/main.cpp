// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include <irrpoly.h>

#include <iostream>
#include <algorithm>

using namespace irrpoly;

[[nodiscard]]
auto generate_irreducible(uintmax_t num) -> std::vector<gfpoly> {
    auto gf2 = make_gf(2);
    // возвращаемое значение
    std::vector<gfpoly> res;
    if (num == 0) {
        return res;
    }
    res.reserve(num);
    res.emplace_back(gfpoly(gf2, {0, 1}));
    if (num == 1) {
        return res;
    }

    // создаём всё необходимое для многопоточности
    multithread::polychecker ch;

    // функция, генерирующая многочлены для проверки
    auto input = [&]() -> gfpoly {
        static uintmax_t index = 1;
        uint8_t degree = 0;
        for (uint8_t i = 1; index >> i; ++i, ++degree);
        std::vector<uintmax_t> res;
        res.reserve(degree + 2);
        res.emplace_back(1);
        for (uint8_t i = 0; i <= degree; ++i) {
            res.emplace_back((index & (1ull << i)) ? 1 : 0);
        }
        ++index;
        return gfpoly(gf2, res);
    };

    auto check = multithread::make_check_func(
        multithread::irreducible_method::berlekamp,
        multithread::primitive_method::nil);

    uintmax_t n = num - 1;
    // функция, вызываемая по окончании проверки, если результат нам подходит - сохраняем и возвращаем true, иначе false
    auto callback = [&](const gfpoly &poly, const typename multithread::result_type &result) -> bool {
        if (result.irreducible) {
            --n;
            res.emplace_back(poly);
        }
        return !n;
    };

    // последний false говорит, что нам нужны все результаты проверки, включая лишние, поскольку мы загружаем
    // многочлены на проверку последовательно, но многопоточность не гарантирует строгого порядка, и какие-то
    // многочлены из начала последовательности могли провериться только после того, как мы уже набрали необходимое
    // количество подходящих, поэтому нужно обработать и их, т.е. последовательность будет избыточна
    ch.pipe(input, check, callback, false);

    // сортируем многочлены в лексико-графическом порядке для получения правильной последовательности
    std::sort(res.begin(), res.end(), [](const gfpoly &a, const gfpoly &b) {
        if (a.degree() == b.degree()) {
            for (auto i = a.degree(); i > 0; --i) {
                if (a[i] != b[i]) {
                    return a[i] < b[i];
                }
            }
            return a[0] < b[0];
        }
        return a.degree() < b.degree();
    });
    // выкидываем лишние с конца, чтобы осталось только требуемое число
    while (res.size() > num) {
        res.pop_back();
    }

    return res;
}

auto main() -> int {
    // генерируем многочлены и выводим результат
    auto poly = generate_irreducible(5);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }
    return 0;
}
