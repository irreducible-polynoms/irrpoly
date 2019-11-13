#include <iostream>
#include <thread>

#include "checker.hpp"

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
        const typename checker<P>::irreducible_method irr_meth,
        const typename checker<P>::primitive_method prim_meth,
        const unsigned threads_num
) {
    // возвращаемое значение
    std::vector<polynomialgf<P>> arr;
    arr.reserve(num);

    // создаём всё необходимое для многопоточности
    checker<P> ch(irr_meth, prim_meth, threads_num);

    auto input = [&]() -> polynomialgf<P> {
        return random<P>(degree);
    };

    auto callback = [&](const typename checker<P>::result_type &res, const polynomialgf<P> &poly) -> bool {
        if (res.irreducible) {
            arr.emplace_back(poly);
            return true;
        }
        return false;
    };

    ch.check(num, input, callback);

    return arr;
}

/// Функция, формирующая данные для тестов. Результаты скопировать в файл test_input.
void test() {
    {
        const uint32_t P = 2;
        for (size_t i = 4; i < 24; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random<P>(i);
                std::cout << (is_irreducible_rabin(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 3;
        for (size_t i = 4; i < 18; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random<P>(i);
                std::cout << (is_irreducible_rabin(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 5;
        for (size_t i = 4; i < 14; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random<P>(i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 7;
        for (size_t i = 4; i < 10; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random<P>(i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 11;
        for (size_t i = 4; i < 6; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random<P>(i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 13;
        for (size_t i = 4; i < 4; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random<P>(i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
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
    const typename checker<P>::irreducible_method irr_meth = checker<P>::irreducible_method::berlekamp;
    // какой из методов проверки на примитивность хотим использовать
    // возможные варианты - отсутствие приверки (nil) и проверка по определению (definition)
    const typename checker<P>::primitive_method prim_meth = checker<P>::primitive_method::nil;
    // число потоков, выполняюих проверку многочленов на соответствие заданной характеристике
    // должно быть равно числу физических ядер минус один
    const unsigned threads_num = std::thread::hardware_concurrency() - 1;

    // генерируем многочлены и выводим результат
    auto poly = generate_irreducible<P>(num, degree, irr_meth, prim_meth, threads_num);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }

//    std::cout << "TEST DATA" << std::endl;
//    test();
    // copy output and save to file named "test_input" (WITHOUT EXTENSION!)
    // use Wolfram Mathematica (program inside test_prog.nb)
    // define function for validating the result
    // f[irr_, prim_, m_, p_] :=
    // irr == Boole[
    //    IrreduciblePolynomialQ[Dot[Power[x, Range[0, Length[p] - 1]], p],
    //     Modulus -> m]] &&
    //  prim == Boole[
    //    PrimitivePolynomialQ[Dot[Power[x, Range[0, Length[p] - 1]], p],
    //     Modulus -> m]]
    // open "input" file for read (CHANGE PATH TO FILE!)
    // file = OpenRead["/Users/vadimpiven/Downloads/test_input"]
    // validate results line by line
    // While[Not[(irr = Read[file, Number]) === EndOfFile],
    // prim = Read[file, Number]; m = Read[file, Number];
    // p = Read[StringToStream[#], Number] & /@
    //   StringSplit[StringTrim[ReadLine[file], ("{" | "}" | " ") ...],
    //    ", "]; Print[f[irr, prim, m, p]]]
    // EXPECTED RESULT IS COLUMN OF "True", if there is some False - this program has some errors
    // close the file
    // Close[file]

    return 0;
}
