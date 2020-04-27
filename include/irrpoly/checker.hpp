/**
 * @file    checker.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include <thread>
#include <cassert>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <utility>
#include <vector>
#include <memory>
#include <optional>
#include <algorithm>

namespace irrpoly {

/// Выполняет проверку на неприводимось и примитивность заданного многочлена над полем GF[P].
template<typename value_type, typename result_type>
class checker {
public:
    /// Вид функции, генерирующей многочлены для проверки.
    using input_func = std::function<value_type()>;

    /// Вид функции, выполняющей проверку и сохраняющей результат.
    using check_func = std::function<void(const value_type &, std::optional<result_type> &)>;

    /// Вид функции, обрабатывающей результат проверки (если многочлен удовлетворяет
    /// требуемым условиям возвращает true, иначе - false).
    using callback_func = std::function<bool(const value_type &, const result_type &)>;

private:
    /// Потоки, непосредственно выполняющие проверку многочленов на неприводимость и примитивность.
    class node {
        std::shared_ptr<std::mutex> s_mutex;
        std::shared_ptr<std::condition_variable> s_cond;

        check_func m_cf; ///< Основная функция проверки

        volatile bool m_terminate; ///< Поток должен быть завершён
        std::optional<value_type> m_val; ///< Входные данные

        bool m_busy; ///< Поток занят полезной работой
        std::optional<result_type> m_res; ///< Результат проверки

        std::thread m_thread;
        std::mutex m_mutex;
        std::condition_variable m_cond;

        /**
         * Собственно функция, выполняющая проверку многочлена на неприводимость и примитивность.
         * Поток бесконечно ожидает получения новых данных. Если данные получены - выполняется их
         * проверка, выставление результата и уведомление условной переменной.
         * Кроме того, постоянно проверяется, не должен ли поток завершить работу.
         * Выход из функции прекращает работу потока и освобожает его ресурсы.
         * В случае с pthread для этого требуется выполнить pthread_exit.
         */
        void check() {
            std::unique_lock<std::mutex> lk(m_mutex);

            while (!m_terminate) {
                if (!(m_busy || m_terminate)) {
                    m_cond.wait(lk);
                    continue;
                }

                m_cf(m_val.value(), m_res);

                std::lock_guard<std::mutex> lg(*s_mutex);
                m_busy = false;
                s_cond->notify_one();
            }
        }

    public:
        node(std::shared_ptr<std::mutex> s_mutex, std::shared_ptr<std::condition_variable> s_cond) :
            s_mutex(std::move(s_mutex)), s_cond(std::move(s_cond)), m_busy(false), m_terminate(false),
            m_val(), m_res() {
            m_thread = std::thread(&node::check, std::ref(*this));
        }

        ~node() {
            m_thread.join();
        }

        /// Устанавливает функцию, которая будет вызываться в процессе проверки
        void set_check(check_func cf) {
            m_cf = cf;
        }

        /// Многочлен, проверка которого выполнялась.
        [[nodiscard]]
        auto get() const -> const std::optional<value_type> & {
            return m_val;
        }

        /**
         * Установка нового многочлена для проверки, сбрасывает результаты
         * предыдущей проверки и выставляет busy = true.
         */
        void set(value_type v) {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_val.emplace(std::move(v));
            m_res.reset();
            m_busy = true;
            m_cond.notify_one();
        }

        /// Возвращает текущее состояние потока.
        [[nodiscard]]
        auto busy() const -> bool {
            return m_busy;
        }

        /// Устанавливает флаг, требующий завершить работу потока по завершении вычислений.
        void terminate() {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_terminate = true;
            m_cond.notify_one();
        }

        /// Возвращает результат проверки текущего многочлена на неприводимость и примитивность.
        [[nodiscard]]
        auto result() const -> const std::optional<result_type> & {
            return m_res;
        }
    }; // class node

    std::shared_ptr<std::mutex> s_mutex;
    std::shared_ptr<std::condition_variable> s_cond;

    std::vector<std::unique_ptr<node>> m_nodes;

    /// Считает, сколько потоков заняты.
    auto countBusy() -> unsigned {
        return std::count_if(m_nodes.begin(), m_nodes.end(), std::mem_fn(&node::busy));
    }

public:
    explicit
    checker(unsigned n = std::thread::hardware_concurrency()) {
        s_mutex = std::make_shared<std::mutex>();
        s_cond = std::make_shared<std::condition_variable>();

        if (n > 1) {
            m_nodes.reserve(--n);
            for (unsigned i = 0; i < n; ++i) {
                m_nodes.push_back(std::make_unique<node>(s_mutex, s_cond));
            }
        }
    }

    /// Основной цикл разделения работы на потоки.
    void check(uint64_t n, input_func in, check_func cf, callback_func back, const bool strict = true) {
        if (m_nodes.empty()) {
            while (n > 0) {
                auto input = in();
                std::optional<result_type> result;
                cf(input, result);
                if (back(input, result.value())) {
                    --n;
                }
            }
            return;
        }

        std::unique_lock<std::mutex> lk(*s_mutex);
        std::this_thread::yield();

        // заряжаем многочлены на проверку
        for (const auto &sl : m_nodes) {
            sl->set_check(cf);
            sl->set(in());
        }
        while (n) {
            // ждём свободный поток
            while (countBusy() == m_nodes.size()) {
                s_cond->wait(lk);
            }
            // находим свободные потоки и заряжаем новыми входными данными
            for (unsigned i = 0; i < m_nodes.size() && n; ++i) {
                if (m_nodes[i]->busy()) {
                    continue;
                }
                if (back(m_nodes[i]->get().value(), m_nodes[i]->result().value())) {
                    --n;
                }
                m_nodes[i]->set(in());
            }
        }
        // ожидаем завершения всех потоков
        while (countBusy()) {
            s_cond->wait(lk);
        }
        if (!strict) {
            // обрабатываем все проверенные многочлены, даже если их больше, чем требовалось найти
            for (const auto &sl : m_nodes) {
                back(sl->get().value(), sl->result().value());
            }
        }
    }

    /// Завершаем работу всех потоков.
    ~checker() {
        for (const auto &sl : m_nodes) {
            sl->terminate();
        }
    }
};

} // namespace irrpoly
