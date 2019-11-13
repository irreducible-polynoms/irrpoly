/**
 * @file    checker.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#ifndef IRRPOLY_CHECKER_HPP
#define IRRPOLY_CHECKER_HPP

#include <thread>
#include <cassert>
#include <functional>

#ifdef PTHREAD

#include <pthread.h>

#else

#include <mutex>
#include <condition_variable>

#endif

#include "polynomialgf.hpp"

namespace irrpoly {

    namespace detail {
        /// Класс для управления мьютексом и условной переменной.
        class sync {
        private:
#ifdef PTHREAD
            pthread_mutex_t mutex;
            pthread_cond_t cond;
#else
            ::std::mutex mutex;
            ::std::unique_lock<::std::mutex> lk;
            ::std::condition_variable cond;
#endif

        public:
            /// Инициализация и блокировка мьютекса.
            sync() noexcept(false) : cond(), mutex() {
#ifdef PTHREAD
                pthread_mutexattr_t attr;
                assert(!pthread_mutexattr_init(&attr));
                assert(!pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK));
                assert(!pthread_mutex_init(&mutex, &attr));
                assert(!pthread_cond_init(&cond, nullptr));
                assert(!pthread_mutexattr_destroy(&attr));
                pthread_mutex_lock(&mutex);
#else
                lk = ::std::unique_lock<::std::mutex>(mutex);
#endif
            }

            /// Освобождение ресурсов, здесь происходит разблокировка мьютекса.
            virtual
#ifndef PTHREAD
            ~sync() = default;

#else
            ~sync() noexcept(false) {
                pthread_mutex_unlock(&mutex);
                assert(!pthread_mutex_destroy(&mutex));
                assert(!pthread_cond_destroy(&cond));
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
    class checker {
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

        /// Вид функции, генерирующей многочлены для проверки.
        typedef ::std::function<polynomialgf<P>()> input_func;

        /// Вид функции, обрабатывающей результат проверки (если многочлен удовлетворяет требуемым условиям возвращает true,
        /// иначе - false).
        typedef ::std::function<bool(const result_type &, const polynomialgf<P> &)> callback_func;

    private:
        /// Потоки, непосредственно выполняющие проверку многочленов на неприводимость и примитивность.
        class node : public detail::sync {
            detail::sync *m; ///< Основной объект синхронизации

            irreducible_method irr_meth; ///< Используемый метод проверки на неприводимость
            primitive_method prim_meth; ///< Используемый метод проверки на примитивность

            polynomialgf<P> poly; ///< Проверяемый многочлен.
            result_type res; ///< Результат проверки

            bool _busy; ///< Поток занят полезной работой
            bool _terminate; ///< Поток должен быть завершён

#ifdef PTHREAD
            pthread_t thread;
#else
            ::std::thread thread;
#endif

            /**
             * Собственно функция, выполняющая проверку многочлена на неприводимость и примитивность.
             * Поток бесконечно ожидает получения новых данных. Если данные получены - выполняется их
             * проверка, выставление результата и уведомление условной переменной.
             * Кроме того, постоянно проверяется, не должен ли поток завершить работу.
             * Выход из функции прекращает работу потока и освобожает его ресурсы.
             * В случае с pthread для этого требуется выполнить pthread_exit.
             */
#ifdef PTHREAD

            static
            void *check(void *arg) noexcept {
                auto *sl = static_cast<node *>(arg);
#else
                void check() noexcept {
                    auto *sl = this;
#endif
                while (!sl->_terminate) {
                    while (!sl->_busy && !sl->_terminate) { sl->wait(); }
                    if (sl->_terminate) { break; }

                    // в случае, когда проверка не выполняется устанавливается результат true
                    switch (sl->irr_meth) {
                        case irreducible_method::berlekamp:
                            sl->res.irreducible = is_irreducible_berlekamp(sl->poly);
                            break;
                        case irreducible_method::rabin:
                            sl->res.irreducible = is_irreducible_rabin(sl->poly);
                            break;
                        default: // irreducible_method::nil
                            sl->res.irreducible = true;
                            break;
                    }
                    switch (sl->prim_meth) {
                        case primitive_method::definition:
                            sl->res.primitive = sl->res.irreducible ? is_primitive_definition(sl->poly) : false;
                            break;
                        default: // primitive_method::nil
                            sl->res.primitive = true;
                            break;
                    }

                    sl->m->lock();
                    sl->_busy = false;
                    sl->m->signal();
                    sl->m->unlock();
                }
#ifdef PTHREAD
                pthread_exit(nullptr);
#endif
            }

        public:
            node(detail::sync *m, irreducible_method irr_meth, primitive_method prim_meth) noexcept :
                    m(m), irr_meth(irr_meth), prim_meth(prim_meth), poly(), res(), _busy(false), _terminate(false),
                    thread() {
#ifdef PTHREAD
                pthread_attr_t thread_attr;
                assert(!pthread_attr_init(&thread_attr));
                assert(!pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED));
                pthread_create(&thread, &thread_attr, &node::check, this);
                assert(!pthread_attr_destroy(&thread_attr));
#else
                thread = ::std::thread(&node::check, this);
                thread.detach();
#endif
            }

            /// Многочлен, проверка которого выполнялась.
            [[nodiscard]]
            const polynomialgf<P> &get() const {
                return poly;
            }

            /**
             * Установка нового многочлена для проверки, сбрасывает результаты
             * предыдущей проверки и выставляет busy = true.
             */
            void set(polynomialgf<P> val) {
                lock();
                poly = val;
                _busy = true;
                signal();
                unlock();
            }

            /// Возвращает текущее состояние потока.
            [[nodiscard]]
            bool busy() const noexcept {
                return _busy;
            }

            /// Устанавливает флаг, требующий завершить работу потока по завершении вычислений.
            void terminate() noexcept {
                lock();
                _terminate = true;
                signal();
                unlock();
            }

            /// Возвращает результат проверки текущего многочлена на неприводимость и примитивность.
            [[nodiscard]]
            const typename checker<P>::result_type &result() const noexcept {
                return res;
            }
        };

        detail::sync m; // master
        ::std::vector<node *> s; // slave

        /// Считает, сколько потоков заняты.
        unsigned countBusy() {
            unsigned res = 0;
            for (const auto *n : s) { res += n->busy(); }
            return res;
        }

    public:
        checker(
                const irreducible_method irr_meth,
                const primitive_method prim_meth,
                const unsigned n = ::std::thread::hardware_concurrency() - 1
        ) noexcept : m(), s() {
            s.reserve(n);
            for (unsigned i = 0; i < n; ++i) {
                s.emplace_back(new node(&m, irr_meth, prim_meth));
            }
        }

        /// Основной цикл разделения работы на потоки.
        void check(uint64_t n, input_func input, callback_func callback, const bool strict = true) noexcept {
            // заряжаем многочлены на проверку
            for (auto *sl : s) { sl->set(input()); }
            while (n) {
                // ждём свободный поток
                while (countBusy() == s.size()) { m.wait(); }
                // находим свободные потоки и заряжаем новыми входными данными
                for (unsigned i = 0; i < s.size() && n; ++i) {
                    if (s[i]->busy()) { continue; }
                    if (callback(s[i]->result(), s[i]->get())) { --n; }
                    s[i]->set(input());
                }
            }
            // ожидаем завершения всех потоков
            while (countBusy()) { m.wait(); }
            if (!strict) {
                // обрабатываем все проверенные многочлены, даже если их больше, чем требовалось найти
                for (auto *sl : s) { callback(sl->result(), sl->get()); }
            }
        }

        /// Завершаем работу всех потоков.
        ~checker() noexcept {
            for (auto *sl : s) { sl->terminate(); }
        }
    };

}

#endif //IRRPOLY_CHECKER_HPP
