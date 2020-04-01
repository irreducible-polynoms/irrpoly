#include <iostream>
#include <thread>
#include <algorithm>

#include <irrpoly/gfcheck.hpp>

using namespace irrpoly;

[[nodiscard]]
std::vector<gfpoly> generate_irreducible(uint64_t num) {
    auto gf2 = make_gf(2);
    // возвращаемое значение
    std::vector<gfpoly> res;
    if (num == 0) { return res; }
    res.reserve(num);
    res.emplace_back(gfpoly(gf2, {0, 1}));
    if (num == 1) { return res; }

    // создаём всё необходимое для многопоточности
    multithread::polychecker ch;

    // функция, генерирующая многочлены для проверки
    auto input = [&]() -> gfpoly {
        static uint64_t index = 1;
        gfpoly res(gf2);
        uint8_t degree = 0;
        for (uint8_t i = 1; index >> i; ++i, ++degree);
        res.data().reserve(degree + 2);
        res.data().emplace_back(gf2, 1);
        for (uint8_t i = 0; i <= degree; ++i) {
            res.data().emplace_back(gf2, (index & (1ull << i)) ? 1 : 0);
        }
        ++index;
        return res;
    };

    auto check = multithread::make_check_func(
            multithread::irreducible_method::berlekamp,
            multithread::primitive_method::nil);

    // функция, вызываемая по окончании проверки, если результат нам подходит - сохраняем и возвращаем true, иначе false
    auto callback = [&](const gfpoly &poly, const typename multithread::result_type &result) -> bool {
        if (!result.irreducible) { return false; }
        res.emplace_back(poly);
        return true;
    };

    // запускаем генерацию num - 1 многочленов (т.к. один у нас уже есть)
    // последний false говорит, что нам нужны все результаты проверки, включая лишние, поскольку мы загружаем
    // многочлены на проверку последовательно, но многопоточность не гарантирует строгого порядка, и какие-то
    // многочлены из начала последовательности могли провериться только после того, как мы уже набрали необходимое
    // количество подходящих, поэтому нужно обработать и их, т.е. последовательность будет избыточна
    ch.check(num - 1, input, check, callback, false);

    // сортируем многочлены в лексико-графическом порядке для получения правильной последовательности
    std::sort(res.begin(), res.end(), [](const gfpoly &a, const gfpoly &b) {
        if (a.degree() == b.degree()) {
            for (uint64_t i = a.degree(); i > 0; --i) {
                if (a[i] != b[i]) { return a[i] < b[i]; }
            }
            return a[0] < b[0];
        } else { return a.degree() < b.degree(); }
    });
    // выкидываем лишние с конца, чтобы осталось только требуемое число
    while (res.size() > num) { res.pop_back(); }

    return res;
}

int main() {
    // генерируем многочлены и выводим результат
    auto poly = generate_irreducible(5);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }

    return 0;
}
