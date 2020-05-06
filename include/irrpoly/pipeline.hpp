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

namespace irrpoly::multithread {

/// Выполняет цепочку "получение входных данных - обработка - возврат результата" в несколько потоков.
template<typename value_t, typename result_t>
class pipeline {
public:
    /// Вид функции, генерирующей входные данные.
    using input_fn = std::function<value_t()>;

    /// Вид функции, выполняющей основную работу и возвращающей её результат.
    using payload_fn = std::function<void(const value_t &, std::optional<result_t> &)>;

    /// Вид функции, обрабатывающей результат и возвращающей true, если работу необходимо завершить.
    using callback_fn = std::function<bool(const value_t &, const result_t &)>;

private:
    /// Потоки, непосредственно выполняющие проверку многочленов на неприводимость и примитивность.
    class pod {
        std::shared_ptr<std::mutex> s_mutex;
        std::shared_ptr<std::condition_variable> s_cond;

        payload_fn m_cf; ///< Основная функция проверки

        volatile bool m_terminate; ///< Поток должен быть завершён
        std::optional<value_t> m_val; ///< Входные данные

        bool m_busy; ///< Поток занят полезной работой
        std::optional<result_t> m_res; ///< Результат проверки

        std::thread m_thread;
        std::mutex m_mutex;
        std::condition_variable m_cond;

        /**
         * Собственно функция, выполняющая основную работу. Поток бесконечно
         * ожидает получения новых данных. Если данные получены - выполняется их
         * проверка, выставление результата и уведомление условной переменной.
         * Кроме того, постоянно проверяется, не должен ли поток завершить работу.
         * Выход из функции прекращает работу потока и освобожает его ресурсы.
         * В случае с pthread для этого требуется выполнить pthread_exit.
         */
        void execute() {
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
        pod(std::shared_ptr<std::mutex> s_mutex, std::shared_ptr<std::condition_variable> s_cond) :
            s_mutex(std::move(s_mutex)), s_cond(std::move(s_cond)), m_busy(false), m_terminate(false),
            m_val(), m_res() {
            m_thread = std::thread(&pod::execute, std::ref(*this));
        }

        ~pod() {
            // Необходимо дождаться завершения потока и лишь затем освобождать
            // память объекта, иначе могут вознивать ошибки доступа к уже
            // освобождённой памяти.
            m_thread.join();
        }

        /// Устанавливает функцию, которая выполняет основную работу
        void set_payload(payload_fn cf) {
            m_cf = cf;
        }

        /// Многочлен, проверка которого выполнялась.
        [[nodiscard]]
        auto get_data() const -> const std::optional<value_t> & {
            return m_val;
        }

        /**
         * Установка новых входных данных, сбрасывает предыдущий результат
         * и выставляет busy = true.
         */
        void set_data(value_t v) {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_val.emplace(std::move(v));
            m_res.reset();
            m_busy = true;
            m_cond.notify_one();
        }

        /// Возвращает текущее состояние потока.
        [[nodiscard]]
        auto is_busy() const -> bool {
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
        auto get_result() const -> const std::optional<result_t> & {
            return m_res;
        }

        /// Очищает хранимые в объекте данные.
        void clear() {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_val.reset();
            m_res.reset();
            m_busy = false;
        }
    }; // class node

    std::shared_ptr<std::mutex> s_mutex;
    std::shared_ptr<std::condition_variable> s_cond;

    std::vector<std::unique_ptr<pod>> m_pods;

    /// Считает, сколько потоков заняты.
    auto countBusy() -> unsigned {
        return std::count_if(m_pods.begin(), m_pods.end(), std::mem_fn(&pod::is_busy));
    }

public:
    explicit
    pipeline(unsigned n = std::thread::hardware_concurrency()) {
        s_mutex = std::make_shared<std::mutex>();
        s_cond = std::make_shared<std::condition_variable>();

        if (n > 1) {
            m_pods.reserve(--n);
            for (unsigned i = 0; i < n; ++i) {
                m_pods.push_back(std::make_unique<pod>(s_mutex, s_cond));
            }
        }
    }

    /// Основной цикл разделения работы на потоки.
    void pipe(input_fn in, payload_fn pl, callback_fn bk, const bool strict = true) {
        if (m_pods.empty()) {
            while (true) {
                auto input = in();
                std::optional<result_t> result;
                pl(input, result);
                if (bk(input, result.value())) {
                    return;
                }
            }
        }

        std::unique_lock<std::mutex> lk(*s_mutex);
        std::this_thread::yield();

        // заряжаем многочлены на проверку
        for (const auto &sl : m_pods) {
            sl->set_payload(pl);
            sl->set_data(in());
        }
        while (true) {
            // ждём свободный поток
            while (countBusy() == m_pods.size()) {
                s_cond->wait(lk);
            }
            // находим свободные потоки и заряжаем новыми входными данными
            for (unsigned i = 0; i < m_pods.size(); ++i) {
                if (m_pods[i]->is_busy()) {
                    continue;
                }
                if (bk(m_pods[i]->get_data().value(), m_pods[i]->get_result().value())) {
                    m_pods[i]->clear();
                    goto end;
                }
                m_pods[i]->set_data(in());
            }
        }
        end:
        // ожидаем завершения всех потоков
        while (countBusy()) {
            s_cond->wait(lk);
        }
        if (!strict) {
            // обрабатываем все проверенные многочлены, даже если их больше, чем требовалось найти
            for (const auto &sl : m_pods) {
                if (sl->get_data().has_value() && sl->get_result().has_value()) {
                    bk(sl->get_data().value(), sl->get_result().value());
                    sl->clear();
                }
            }
        }
    }

    /// Завершаем работу всех потоков.
    ~pipeline() {
        for (const auto &sl : m_pods) {
            sl->terminate();
        }
    }
};

} // namespace irrpoly::multithread
