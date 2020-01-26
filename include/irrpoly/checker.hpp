/**
 * @file    checker.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#ifndef IRRPOLY_CHECKER_HPP
#define IRRPOLY_CHECKER_HPP

#include <iostream>

#include <thread>
#include <cassert>
#include <functional>
#include <memory>
#include <utility>
#include <list>

#ifdef PTHREAD

#include <pthread.h>

#else

#include <shared_mutex>
#include <condition_variable>

#endif

namespace irrpoly {

    namespace detail {

        /// Класс для управления мьютексом и условной переменной.
        class sync {
        private:
#ifdef PTHREAD
            pthread_mutex_t mutex;
            pthread_cond_t cond;
#else
            ::std::shared_mutex mutex;
            ::std::unique_lock<::std::shared_mutex> lk;
            ::std::condition_variable_any cond;
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
                lk = ::std::unique_lock<::std::shared_mutex>(mutex);
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
        }; // class sync

        /// Класс, обеспечивающий многопоточный доступ к данным
        template<typename T>
        class cell {
        private:
            sync s; ///< Объект синхронизации
            T data; ///< Охраняемые данные
            T def; ///< Значение data по умолчанию, устанавливаемое функцией reset
            bool n; ///< True, если данные были обновлены но ещё не прочитаны

        public:
            /// Вид функции, проверяющей что ожидаемое изменение данных произошло
            typedef ::std::function<bool(const T &)> condition_func;

            /// Владение объектом data передаётся внутрь cell
            explicit
            cell(T data, T def = T()) noexcept : s(), data(std::move(data)), def(std::move(def)), n(false) {}

            /// Возвращает замороженный объект
            [[nodiscard]]
            const T & get() noexcept {
                n = false;
                return data;
            }

            /// Сигнализирует, что есть новые непрочитанные данные
            [[nodiscard]]
            bool updated() const noexcept {
                return n;
            }

            /// Устанавливает значение по умолчанию
            void reset() noexcept {
                n = true;
                data = T(def);
            }

            /// Ждёт разрешения на изменение объекта, меняет его и уведомляет ожидающего
            void set(T d) noexcept {
                s.lock();
                data = std::move(d);
                n = true;
                s.signal();
                s.unlock();
            }

            /// Даёт разрешение на изменение объекта и ждёт, пока изменение произойдёт
            const T & wait_change(condition_func cond) noexcept {
                while (!(n && cond(data))) {
                    s.wait();
                }
                return data;
            }
        }; // class cell

        /// Сетка из ячеек, защищённых общим мьютексом
        template<typename T>
        class grid {
        public:
            class cell {
            private:
                ::std::shared_ptr<sync> s; ///< Объект синхронизации
                T data; ///< Охраняемые данные
                T def; ///< Значение data по умолчанию, устанавливаемое функцией reset

            public:
                /// Владение объектом data передаётся внутрь cell
                explicit
                cell(::std::shared_ptr<sync> s, T data, T def) noexcept :
                    s(std::move(s)), data(std::move(data)), def(std::move(def)) {}

                /// Возвращает замороженный объект
                [[nodiscard]]
                const T & get() const noexcept {
                    return data;
                }

                /// Устанавливает значение по умолчанию
                void reset() noexcept {
                    data = T(def);
                }

                /// Ждёт разрешения на изменение объекта, меняет его и уведомляет ожидающего
                void set(T d) noexcept {
                    s->lock();
                    data = std::move(d);
                    s->signal();
                    s->unlock();
                }
            }; // class grid::cell

        private:
            ::std::shared_ptr<sync> s;
            ::std::list<grid::cell> data;

        public:
            /// Вид функции, проверяющей что ожидаемое изменение данных произошло
            typedef ::std::function<bool(const ::std::list<grid::cell> &)> condition_func;

            grid() noexcept : data() {
                s = ::std::make_shared<sync>();
            }

            /// Добавляем в сетку ещё один элемент
            grid::cell & chain(T d, T def = T()) noexcept {
                return data.emplace_back(s, std::move(d), std::move(def));
            }

            /// Даёт разрешение на изменение объекта и ждёт, пока изменение произойдёт
            void wait_change(condition_func cond) noexcept {
                while (!cond(data)) {
                    s->wait();
                }
            }
        }; // class grid

    } // namespace detail

    /// Выполняет проверку на неприводимось и примитивность заданного многочлена над полем GF[P].
    template<typename value_type, typename result_type>
    class checker {
    public:
        /// Вид функции, генерирующей многочлены для проверки.
        typedef ::std::function<value_type()> input_func;

        /// Вид функции, выполняющей проверку и сохраняющей результат.
        typedef ::std::function<void(const value_type &, result_type &)> check_func;

        /// Вид функции, обрабатывающей результат проверки (если многочлен удовлетворяет
        /// требуемым условиям возвращает true, иначе - false).
        typedef ::std::function<bool(const value_type &, const result_type &)> callback_func;

    private:
        /// Потоки, непосредственно выполняющие проверку многочленов на неприводимость и примитивность.
        class node {
        private:
            typedef struct {
                bool terminate; ///< Поток должен быть завершён
                value_type val; ///< Входные данные
            } outside_ctrl;

            // Управляется снаружи
            detail::cell<outside_ctrl> input;

            // Управляется изнутри
            typename detail::grid<bool>::cell &_busy; ///< Поток занят полезной работой
            result_type res; ///< Результат проверки

            // Инициализируется при конструировании
            check_func cf; ///< Основная функция проверки
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
                while (!sl->input.get().terminate) {
                    sl->input.wait_change([&](const outside_ctrl &data) {
                       return _busy.get() || data.terminate;
                    });
                    if (sl->input.get().terminate) { break; }

                    sl->cf(sl->input.get().val, sl->res);
                    _busy.set(false);
                }
#ifdef PTHREAD
                pthread_exit(nullptr);
#endif
            }

        public:
            explicit
            node(typename detail::grid<bool>::cell &busy) noexcept :
                input({ false, value_type() }), _busy(busy), res(), cf(), thread() {
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

            /// Устанавливает функцию, которая будет вызываться в процессе проверки
            void set_check(check_func c) noexcept {
                cf = c;
            }

            /// Многочлен, проверка которого выполнялась.
            [[nodiscard]]
            const value_type &get() {
                return input.get().val;
            }

            /**
             * Установка нового многочлена для проверки, сбрасывает результаты
             * предыдущей проверки и выставляет busy = true.
             */
            void set(value_type v) {
                _busy.reset();
                input.set({ false, std::move(v) });
            }

            /// Возвращает текущее состояние потока.
            [[nodiscard]]
            bool busy() const noexcept {
                return _busy.get();
            }

            /// Устанавливает флаг, требующий завершить работу потока по завершении вычислений.
            void terminate() noexcept {
                input.set({ true, value_type() });
            }

            /// Возвращает результат проверки текущего многочлена на неприводимость и примитивность.
            [[nodiscard]]
            const result_type &result() const noexcept {
                return res;
            }
        }; // class node

        detail::grid<bool> m; // master
        ::std::vector<::std::unique_ptr<node>> s; // slave

        /// Считает, сколько потоков заняты.
        unsigned countBusy() {
            unsigned res = 0;
            for (const auto &n : s) { res += n.busy(); }
            return res;
        }

    public:
        explicit
        checker(const unsigned n = ::std::thread::hardware_concurrency() - 1) noexcept : m(), s() {
            s.reserve(n);
            for (unsigned i = 0; i < n; ++i) {
                s.push_back(::std::make_unique<node>(m.chain(false, true)));
            }
        }

        /// Основной цикл разделения работы на потоки.
        void check(uint64_t n, input_func in, check_func cf, callback_func back, const bool strict = true) noexcept {
            // заряжаем многочлены на проверку
            for (const auto &sl : s) {
                sl->set_check(cf);
                sl->set(in());
            }
            while (n) {
                // ждём свободный поток
                m.wait_change([&](const ::std::list<detail::grid<bool>::cell> &c) -> bool {
                    unsigned res = 0;
                    for (const auto &n : c) { res += n.get(); }
                    return res < s.size();
                });
                // находим свободные потоки и заряжаем новыми входными данными
                for (unsigned i = 0; i < s.size() && n; ++i) {
                    if (s[i]->busy()) { continue; }
                    if (back(s[i]->get(), s[i]->result())) { --n; }
                    s[i]->set(in());
                }
            }
            // ожидаем завершения всех потоков
            m.wait_change([](const ::std::list<detail::grid<bool>::cell> &c) -> bool {
                for (const auto &n : c) {
                    if (n.get()) return false;
                }
                return true;
            });
            if (!strict) {
                // обрабатываем все проверенные многочлены, даже если их больше, чем требовалось найти
                for (const auto &sl : s) {
                    back(sl->get(), sl->result());
                }
            }
        }

        /// Завершаем работу всех потоков.
        ~checker() noexcept {
            for (const auto &sl : s) {
                sl->terminate();
            }
        }
    }; // class checker

} // namespace irrpoly

#endif //IRRPOLY_CHECKER_HPP
