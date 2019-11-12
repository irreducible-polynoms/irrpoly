/**
 * @file    checker.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#ifndef IRRPOLY_CHECKER_HPP
#define IRRPOLY_CHECKER_HPP

#ifdef PTHREAD

#include <pthread.h>

#else
#include <thread>
#include <mutex>
#include <condition_variable>
#endif

#include "polynomialgf.hpp"

namespace detail {
    /// Класс для управления мьютексом и условной переменной.
    class sync {
    private:
#ifdef PTHREAD
        pthread_mutex_t mutex;
        pthread_cond_t cond;
#else
        std::mutex mutex;
        std::unique_lock<std::mutex> lk;
        std::condition_variable cond;
#endif

    public:
        /// Инициализация и блокировка мьютекса.
        sync() noexcept(false) : cond(), mutex() {
#ifdef PTHREAD
            pthread_mutexattr_t attr;
            if (pthread_mutexattr_init(&attr) ||
                pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK) ||
                pthread_mutex_init(&mutex, &attr) ||
                pthread_cond_init(&cond, nullptr) ||
                pthread_mutexattr_destroy(&attr)) {
                throw std::runtime_error("pthread init failed");
            }
            pthread_mutex_lock(&mutex);
#else
            lk = std::unique_lock<std::mutex>(mutex);
#endif
        }

#ifdef PTHREAD

        /// Освобождение ресурсов, здесь происходит разблокировка мьютекса.
        virtual
        ~sync() noexcept(false) {
            pthread_mutex_unlock(&mutex);
            if (pthread_mutex_destroy(&mutex) ||
                pthread_cond_destroy(&cond)) {
                throw std::runtime_error("pthread destroy failed");
            }
        }

#endif

        /// Блокировка мьютекса.
        void lock() noexcept {
#ifdef PTHREAD
            pthread_mutex_lock(&mutex);
#else
            mutex.lock();
#endif
        }

        /// Разблокировка мьютекса.
        void unlock() noexcept {
#ifdef PTHREAD
            pthread_mutex_unlock(&mutex);
#else
            mutex.unlock();
#endif
        }

        /// Ожидание уведомления условной переменной.
        void wait() noexcept {
#ifdef PTHREAD
            pthread_cond_wait(&cond, &mutex);
#else
            cond.wait(lk);
#endif
        }

        /// Уведомление условной переменной.
        void signal() noexcept {
#ifdef PTHREAD
            pthread_cond_signal(&cond);
#else
            cond.notify_one();
#endif
        }
    };
}


/// Выполняет проверку на неприводимось и примитивность заданного многочлена над полем GF[P].
template<uint32_t P>
class checker : private detail::sync {
public:
    /// Структура, представляющая результаты проверки многочлена.
    struct result_type {
        bool irreducible;
        bool primitive;
    };

    /// Доступные методы проверки нанеприводимость.
    enum class irreducible_method {
        nil, ///< не проверять
        berlekamp, ///< алгоритм Берлекампа
        rabin ///< алгоритм Рабина
    };

    /// Доступные методы проверки примитивность.
    enum class primitive_method {
        nil, ///< не проверять
        definition ///< проверка по определению
    };

    /// Cлужебный класс, используется для выполнения проверки в несколько потоков.
    class control_type : public detail::sync {
    private:
        /**
         * Для проверки требуется мьютекс, условная переменная, заданное число потоков
         * и по одному хранилищу данных на каждый поток.
         */
#ifdef PTHREAD
        std::vector<pthread_t> threads;
#else
        std::vector<std::thread> threads;
#endif
        std::vector<checker<P> *> _checkers;

    public:
        /**
         * Инициализация мьютекса и условной переменной, создание требуемого числа потоков
         * и хранилищ данных. На выходе мьютекс находится в состоянии locked, все потоки
         * запущены и находятся в состоянии detached.
         */
        control_type(
                const typename checker<P>::irreducible_method irr_meth,
                const typename checker<P>::primitive_method prim_meth,
                const unsigned threads_num
        ) noexcept(false) {
			_checkers.reserve(threads_num);
			for (unsigned i = 0; i < threads_num; ++i) {
				_checkers.emplace_back(new checker<P>(this, irr_meth, prim_meth));
			}
#ifdef PTHREAD
            threads = std::vector<pthread_t>(threads_num);

            pthread_attr_t thread_attr;
            if (pthread_attr_init(&thread_attr) ||
                pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED)) {
                throw std::runtime_error("pthread init failed");
            }

            for (unsigned i = 0; i < threads_num; ++i) {
                pthread_create(&threads[i], &thread_attr, &checker<P>::check, _checkers[i]);
            }

            if (pthread_attr_destroy(&thread_attr)) {
                throw std::runtime_error("pthread destroy failed");
            }
#else
            threads.reserve(threads_num);
            for (unsigned i = 0; i < threads_num; ++i) {
                threads.emplace_back(std::thread(&checker<P>::check, _checkers[i]));
                threads.back().detach();
            }
#endif
        }

        /**
         * Возвращает структуры, хранящие все данные потока: проверяемый многочлен,
         * текущее состояение потока и результат проверки многочлена.
         */
        std::vector<checker<P>*> &checkers() noexcept {
            return _checkers;
        }

        ~control_type() noexcept {
            for (auto *c : _checkers) {
                delete c;
            }
        }
    };

private:
    /// Проверяемый многочлен.
    polynomialgf<P> poly;

    /// Ссылка на управляющую структуру, нужна для использования мьютекса и условной переменной.
    control_type *ctrl;

    bool _busy; ///< поток занят полезной работой
    bool _terminate; ///< поток должен быть завершён
    result_type res; ///< результат проверки
    irreducible_method irr_meth; ///< используемый метод проверки на неприводимость
    primitive_method prim_meth; ///< используемый метод проверки на примитивность

public:
    explicit
    checker(control_type *, irreducible_method, primitive_method) noexcept;

    polynomialgf<P> get() const;

    void set(polynomialgf<P>);

    [[nodiscard]]
    bool busy() const noexcept;

    void terminate() noexcept;

    [[nodiscard]]
    const result_type &result() const noexcept;

#ifdef PTHREAD

    static
    void *check(void *arg) noexcept;

#else
    void check() noexcept;
#endif
};

/**
 * Конструктор структуры с данными для проверки и результатом.
 * Начальное состояние - busy = false.
 */
template<uint32_t P>
checker<P>::checker(
        control_type *ctrl,
        const irreducible_method irr_meth,
        const primitive_method prim_meth
) noexcept :
        poly(), ctrl(ctrl), _busy(false), res({false, false}),
        irr_meth(irr_meth), prim_meth(prim_meth), _terminate(false) {}

/// Многочлен, проверка которого выполнялась.
template<uint32_t P>
polynomialgf<P> checker<P>::get() const {
    return poly;
}

/**
 * Установка нового многочлена для проверки, сбрасывает результаты
 * предыдущей проверки и выставляет busy = true.
 */
template<uint32_t P>
void checker<P>::set(polynomialgf<P> val) {
    lock();
    poly = val;
    _busy = true;
    signal();
    unlock();
}

/**
 * Собственно функция, выполняющая проверку многочлена на неприводимость и примитивность.
 * Поток бесконечно ожидает получения новых данных. Если данные получены - выполняется их
 * проверка, выставление результата и уведомление условной переменной.
 * Кроме того, постоянно проверяется, не должен ли поток завершить работу.
 * Выход из функции прекращает работу потока и освобожает его ресурсы.
 * В случае с pthread для этого требуется выполнить pthread_exit.
 */
template<uint32_t P>
#ifdef PTHREAD
void *checker<P>::check(void *arg) noexcept {
    auto *c = static_cast<checker *>(arg);
#else
    void checker<P>::check() noexcept {
        auto *c = this;
#endif
    while (!c->_terminate) {
        while (!c->_busy && !c->_terminate) { c->wait(); }
        if (c->_terminate) { break; }

        // в случае, когда проверка не выполняется устанавливается результат true
        switch (c->irr_meth) {
            case irreducible_method::berlekamp:
                c->res.irreducible = is_irreducible_berlekamp(c->poly);
                break;
            case irreducible_method::rabin:
                c->res.irreducible = is_irreducible_rabin(c->poly);
                break;
            default: // irreducible_method::nil
                c->res.irreducible = true;
                break;
        }
        switch (c->prim_meth) {
            case primitive_method::definition:
                c->res.primitive = c->res.irreducible ? is_primitive_definition(c->poly) : false;
                break;
            default: // primitive_method::nil
                c->res.primitive = true;
                break;
        }

        c->ctrl->lock();
        c->_busy = false;
        c->ctrl->signal();
        c->ctrl->unlock();
    }
#ifdef PTHREAD
    pthread_exit(nullptr);
#endif
}

/// Возвращает текущее состояние потока.
template<uint32_t P>
[[nodiscard]]
bool checker<P>::busy() const noexcept {
    return _busy;
}

/// Устанавливает флаг, требующий завершить работу потока по завершении вычислений.
template<uint32_t P>
void checker<P>::terminate() noexcept {
    _terminate = true;
    if (!_busy) {
        lock();
        signal();
        unlock();
    }
}

/// Возвращает результат проверки текущего многочлена на неприводимость и примитивность.
template<uint32_t P>
[[nodiscard]]
const typename checker<P>::result_type &checker<P>::result() const noexcept {
    return res;
}

#endif //IRRPOLY_CHECKER_HPP
