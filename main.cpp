#include <iostream>
#include <stdexcept>
#include <thread>

#include <pthread.h>
#include "checker.hpp"
#include "polynomialgf.hpp"

template <uint32_t P>
[[nodiscard]]
polynomialgf<P> generate_primitive(typename polynomialgf<P>::size_type degree, const typename checker<P>::method meth) {
    static const unsigned threadsNum = (std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 1);

    auto countBusy = [](const std::vector<checker<P>> &c) noexcept {
        unsigned n = 0;
        for (auto &checker : c) { n += checker.busy(); }
        return n;
    };

    polynomialgf<P> res; // возвращаемое значение

    // всё для работы с потоками
    pthread_t threads[threadsNum];
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutexattr_t attr;
    pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

    if (pthread_mutexattr_init(&attr) ||
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK) ||
        pthread_mutex_init(&mutex, &attr) ||
        pthread_cond_init(&cond, nullptr)) { throw std::runtime_error("pthread initialisation failed"); }

    // по одному проверщику в каждом потоке
    auto c = std::vector<checker<P>>(threadsNum, checker<P>(&mutex, &cond, meth));

    while (true) {
        while (countBusy(c) >= threadsNum) {
            // ждём изменения числа занятых потоков
            pthread_cond_wait(&cond, &mutex);
        }

        for (unsigned j = 0; j < threadsNum; ++j) {
            // находим первый свободный поток
            if (c[j].busy()) { continue; }
            if (c[j].result().primitive) {
                res = c[j].get();
                goto END;
            }
            // генерируем случайный многочлен для проверки
            c[j].set(random<P>(degree));
            // создаём новый поток для выполнения проверки
            pthread_create(&threads[j], nullptr, &checker<P>::check, &c[j]);
            // отсоединеняем поток
            pthread_detach(threads[j]);
        }
    }
    END:
    // ждём завершения всех созданных потоков
    while (countBusy(c)) {
        pthread_cond_wait(&cond, &mutex);
    }
    if (pthread_cond_destroy(&cond) || pthread_mutex_destroy(&mutex)
        || pthread_mutexattr_destroy(&attr)) { throw std::runtime_error("pthread destruction failed"); }
    return res;
}

int main() {
    const uint32_t P = 2;
    std::cout << generate_primitive<P>(5, checker<P>::method::rabin) << std::endl;

    return 0;
}
