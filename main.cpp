#include <iostream>
#include <thread>

#include "checker.hpp"

/**
 * Шаблон функции, применяемой для генерации многочленов с требуемыми характеристиками.
 * В данном случае генерируется заданное количество примитивных многочленов требуемой степени.
 * В функции генерируются случайные многочлены требуемой степени и выполняется их проверка
 * на неприводимость и примитивность. Для проверки на неприводимость можно использовать
 * один из алгоритмов - Рабина или Берлекампа. Для осмысленного применения рекоммендуется
 * перебирать многочлены последовательно, а не выбирать каждый раз случайный (во избежание
 * повторов). Проверка многочленов выполняется параллельно заданным числом потоков. Число
 * потоков должно соответствовать числу физических ядер минус один (одно ядро задействуется
 * основным потоком приложения, передающим многочлены на проверку в другие потоки).
 */
template<uint32_t P>
[[nodiscard]]
std::vector<polynomialgf<P>> generate_primitive(
        typename std::vector<polynomialgf<P>>::size_type num,
        const typename polynomialgf<P>::size_type degree,
        const typename checker<P>::method meth,
        const unsigned threads_num
) {
    // функция для подсчёта потоков, в данный момент выполняющих проверку
    auto countBusy = [](const std::vector<checker<P>> &c) noexcept {
        unsigned n = 0;
        for (const auto &checker : c) { n += checker.busy(); }
        return n;
    };

    // возвращаемое значение
    std::vector<polynomialgf<P>> res(num--);

    // создаём всё необходимое для многопоточности
    typename checker<P>::control_type ctrl(meth, threads_num);

    // пока не нашли необходимое количество многочленов требуемой характеристики
    while (true) {
        // проверяем, что есть свободный поток
        if (countBusy(ctrl.checkers()) == threads_num) {
            // ждём, пока свободный поток появится
            ctrl.wait();
        }

        // находим свободный поток
        for (unsigned i = 0; i < threads_num; ++i) {
            if (ctrl.checkers()[i].busy()) { continue; }
            // если многочлен примитивный (именно их и ищем) - сохраняем результат
            // для поиска неприводимых используйте .irreducible
            if (ctrl.checkers()[i].result().primitive) {
                res[num] = ctrl.checkers()[i].get();
                if (num > 0) { num -= 1; }
                    // если нашли необходимое число многочленов заданной
                    // характеристики - выходим из цикла
                else { goto END; }
            }
            // генерируем новый многочлен для проверки (в данном случае просто случайный)
            ctrl.checkers()[i].set(random<P>(degree));
        }
    }
    END:
    // сообщаем всем потокам, что они должны прекратить работу
    for (auto &c : ctrl.checkers()) {
        c.terminate();
    }
    // ждём завершения всех потоков
    while (countBusy(ctrl.checkers())) {
        ctrl.wait();
    }

    return res;
}

int main() {
    // основание поля Галуа GF[P]
    const uint32_t P = 2;
    // число многочленов, которые требуется найти
    const typename std::vector<polynomialgf<P>>::size_type num = 3;
    // степень искомых многочленов
    const typename polynomialgf<P>::size_type degree = 5;
    // какой из методов проверки на неприводимость хотим использовать
    // возможные вариенты - Берлекампа (berlekamp) и Рабин (rabin)
    const typename checker<P>::method meth = checker<P>::method::rabin;
    // число потоков, выполняюих проверку многочленов на соответствие заданной характеристике
    // должно быть равно числу физических ядер минус один
    const unsigned threads_num = std::thread::hardware_concurrency() - 1;

    // генерируем многочлены и выводим результат
    auto poly = generate_primitive<P>(num, degree, meth, threads_num);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }

    return 0;
}
