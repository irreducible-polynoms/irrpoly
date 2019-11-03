#include <iostream>
#include <thread>

#include "checker.hpp"

template<uint32_t P>
[[nodiscard]]
std::vector<polynomialgf<P>> generate_primitive(
        typename std::vector<polynomialgf<P>>::size_type num,
        const typename polynomialgf<P>::size_type degree,
        const typename checker<P>::method meth,
        const unsigned threads_num
) {
    auto countBusy = [](const std::vector<checker<P>> &c) noexcept {
        unsigned n = 0;
        for (const auto &checker : c) { n += checker.busy(); }
        return n;
    };

    std::vector<polynomialgf<P>> res(num--); // возвращаемое значение

    // создаём всё необходимое для многопоточности
    typename checker<P>::control_type ctrl;
    checker<P>::init(ctrl, meth, threads_num);
    // закрываем мьютекс, таким образом изменения состояния busy() отдельных потоков
    // на значение false смогут происходить только внутри pthread_cond_wait
    pthread_mutex_lock(&ctrl.mutex);

    while (true) {
        if (countBusy(ctrl.checkers) == threads_num) {
            // ждём изменения числа занятых потоков
            pthread_cond_wait(&ctrl.cond, &ctrl.mutex);
        }

        for (unsigned i = 0; i < threads_num; ++i) {
            // находим первый свободный поток
            if (ctrl.checkers[i].busy()) { continue; }
            // если многочлен примитивный (именно их и ищем) - созраняем результат
            // для поиска неприводимых используйте .irreducible
            if (ctrl.checkers[i].result().primitive) {
                res[num] = ctrl.checkers[i].get();
                if (num > 0) { num -= 1; }
                else { goto END; }
            }
            // генерируем случайный многочлен для проверки
            ctrl.checkers[i].set(random<P>(degree));
        }
    }
    END:
    // сообщаем всем потокам, что они должны прекратить работу
    for (auto &c : ctrl.checkers) {
        c.terminate();
    }
    // ждём завершения всех потоков
    while (countBusy(ctrl.checkers)) {
        pthread_cond_wait(&ctrl.cond, &ctrl.mutex);
    }

    // освобождаем мьютекс, т.к. работа с потоками окончена
    pthread_mutex_unlock(&ctrl.mutex);
    // освобождаем ресурсы
    checker<P>::destroy(ctrl);

    return res;
}

int main() {
    const uint32_t P = 2;
    const typename std::vector<polynomialgf<P>>::size_type num = 3;
    const typename polynomialgf<P>::size_type degree = 5;
    const typename checker<P>::method meth = checker<P>::method::rabin;
    const unsigned threads_num = std::thread::hardware_concurrency() - 1;

    auto poly = generate_primitive<P>(num, degree, meth, threads_num);
    for (const auto &p : poly) {
        std::cout << p << std::endl;
    }

    return 0;
}
