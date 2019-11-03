#include <iostream>
#include <stdexcept>
#include <thread>

#include <pthread.h>
#include "checker.hpp"
#include "polynomialgf.hpp"

template<uint32_t P>
[[nodiscard]]
polynomialgf<P> generate_primitive(
        typename polynomialgf<P>::size_type degree,
        const typename checker<P>::method meth,
        const unsigned threadsNum
) {
    auto countBusy = [](const std::vector<checker<P>> &c) noexcept {
        unsigned n = 0;
        for (auto &checker : c) { n += checker.busy(); }
        return n;
    };

    polynomialgf<P> res; // возвращаемое значение

    // всё для работы с потоками
    pthread_mutex_t mutex;
    pthread_mutexattr_t attr;
    pthread_cond_t cond;

    if (pthread_mutexattr_init(&attr) ||
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK) ||
        pthread_mutex_init(&mutex, &attr) ||
        pthread_cond_init(&cond, nullptr)) {
        throw std::runtime_error("pthread initialisation failed");
    }

    // закрываем мьютекс, таким образом изменения состояния busy() отдельных потоков
    // на значение false смогут происходить только внутри pthread_cond_wait
    pthread_mutex_lock(&mutex);

    // по одному проверщику на каждый поток
    auto checkers = std::vector<checker<P>>(threadsNum, checker<P>(&mutex, &cond, meth));
    std::vector<pthread_t> threads(threadsNum);
    pthread_attr_t thread_attr;
    if (pthread_attr_init(&thread_attr) ||
        pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED)) {
        throw std::runtime_error("pthread initialisation failed");
    }
    for (unsigned i = 0; i < threadsNum; ++i) {
        // создаём новый поток для выполнения проверки
        pthread_create(&threads[i], &thread_attr, &checker<P>::check, &checkers[i]);
    }

    while (true) {
        if (countBusy(checkers) == threadsNum) {
            // ждём изменения числа занятых потоков
            pthread_cond_wait(&cond, &mutex);
        }

        for (unsigned i = 0; i < threadsNum; ++i) {
            // находим первый свободный поток
            if (checkers[i].busy()) { continue; }
            if (checkers[i].result().primitive) {
                res = checkers[i].get();
                goto END;
            }
            // генерируем случайный многочлен для проверки
            checkers[i].set(random<P>(degree));
        }
    }
    END:
    // сообщаем всем потокам, что они должны прекратить работу
    for (auto &c : checkers) {
        c.terminate();
    }
    // ждём завершения всех потоков
    while (countBusy(checkers)) {
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    if (pthread_cond_destroy(&cond) ||
        pthread_mutex_destroy(&mutex) ||
        pthread_mutexattr_destroy(&attr) ||
        pthread_attr_destroy(&thread_attr)) {
        throw std::runtime_error("pthread destruction failed");
    }
    return res;
}

int main() {
    const uint32_t P = 2;
    const unsigned threadsNum = std::thread::hardware_concurrency() - 1;
    std::cout << generate_primitive<P>(5, checker<P>::method::rabin, threadsNum) << std::endl;

    return 0;
}
