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
#include <mutex>
#include <condition_variable>
#include <utility>
#include <vector>
#include <memory>

namespace irrpoly {

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
            ::std::shared_ptr<::std::mutex> s_mutex;
            ::std::shared_ptr<::std::condition_variable> s_cond;

            check_func cf; ///< Основная функция проверки

            bool _terminate; ///< Поток должен быть завершён
            value_type val; ///< Входные данные

            bool _busy; ///< Поток занят полезной работой
            result_type res; ///< Результат проверки

            ::std::thread thread;
            ::std::mutex mutex;
            ::std::condition_variable cond;

            /**
             * Собственно функция, выполняющая проверку многочлена на неприводимость и примитивность.
             * Поток бесконечно ожидает получения новых данных. Если данные получены - выполняется их
             * проверка, выставление результата и уведомление условной переменной.
             * Кроме того, постоянно проверяется, не должен ли поток завершить работу.
             * Выход из функции прекращает работу потока и освобожает его ресурсы.
             * В случае с pthread для этого требуется выполнить pthread_exit.
             */
            void check() {
                ::std::unique_lock<::std::mutex> lk(mutex);

                while (!_terminate) {
                    while (!(_busy || _terminate)) { cond.wait(lk); }
                    if (_terminate) { break; }

                    cf(val, res);

                    ::std::lock_guard<::std::mutex> lg(*s_mutex);
                    _busy = false;
                    s_cond->notify_one();
                }
            }

        public:
            node(::std::shared_ptr<::std::mutex> s_mutex, ::std::shared_ptr<::std::condition_variable> s_cond) :
                    s_mutex(std::move(s_mutex)), s_cond(std::move(s_cond)), _busy(false), _terminate(false) {
                thread = ::std::thread(&node::check, this);
                thread.detach();
            }

            ~node() {
                if (thread.joinable()) { thread.join(); }
            }

            /// Устанавливает функцию, которая будет вызываться в процессе проверки
            void set_check(check_func c) {
                cf = c;
            }

            /// Многочлен, проверка которого выполнялась.
            [[nodiscard]]
            const value_type &get() const {
                return val;
            }

            /**
             * Установка нового многочлена для проверки, сбрасывает результаты
             * предыдущей проверки и выставляет busy = true.
             */
            void set(value_type v) {
                ::std::lock_guard<::std::mutex> lg(mutex);
                val = v;
                _busy = true;
                cond.notify_one();
            }

            /// Возвращает текущее состояние потока.
            [[nodiscard]]
            bool busy() const {
                return _busy;
            }

            /// Устанавливает флаг, требующий завершить работу потока по завершении вычислений.
            void terminate() {
                ::std::lock_guard<::std::mutex> lg(mutex);
                _terminate = true;
                cond.notify_one();
            }

            /// Возвращает результат проверки текущего многочлена на неприводимость и примитивность.
            [[nodiscard]]
            const result_type &result() const {
                return res;
            }
        }; // class node

        ::std::shared_ptr<::std::mutex> s_mutex;
        ::std::shared_ptr<::std::condition_variable> s_cond;

        ::std::vector<::std::unique_ptr<node>> s; // slave

        /// Считает, сколько потоков заняты.
        unsigned countBusy() {
            unsigned res = 0;
            for (const auto &n : s) { res += n->busy(); }
            return res;
        }

    public:
        explicit
        checker(const unsigned n = ::std::thread::hardware_concurrency() - 1) {
            s_mutex = ::std::make_shared<::std::mutex>();
            s_cond = ::std::make_shared<::std::condition_variable>();

            s.reserve(n);
            for (unsigned i = 0; i < n; ++i) {
                s.push_back(::std::make_unique<node>(s_mutex, s_cond));
            }
        }

        /// Основной цикл разделения работы на потоки.
        void check(uint64_t n, input_func in, check_func cf, callback_func back, const bool strict = true) {
            ::std::unique_lock<::std::mutex> lk(*s_mutex);
            // заряжаем многочлены на проверку
            for (const auto &sl : s) {
                sl->set_check(cf);
                sl->set(in());
            }
            while (n) {
                // ждём свободный поток
                while (countBusy() == s.size()) { s_cond->wait(lk); }
                // находим свободные потоки и заряжаем новыми входными данными
                for (unsigned i = 0; i < s.size() && n; ++i) {
                    if (s[i]->busy()) { continue; }
                    if (back(s[i]->get(), s[i]->result())) { --n; }
                    s[i]->set(in());
                }
            }
            // ожидаем завершения всех потоков
            while (countBusy()) { s_cond->wait(lk); }
            if (!strict) {
                // обрабатываем все проверенные многочлены, даже если их больше, чем требовалось найти
                for (const auto &sl : s) { back(sl->get(), sl->result()); }
            }
        }

        /// Завершаем работу всех потоков.
        ~checker() {
            for (const auto &sl : s) { sl->terminate(); }
        }
    };

} // namespace irrpoly

#endif //IRRPOLY_CHECKER_HPP