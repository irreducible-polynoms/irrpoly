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

template<uint32_t P>
class checker {
public:
    struct result_type {
        bool irreducible;
        bool primitive;
    };

    enum method {
        berlekamp = 0,
        rabin
    };

    class control_type {
    private:
#ifdef PTHREAD
        pthread_mutex_t mutex;
        pthread_cond_t cond;
        std::vector<pthread_t> threads;
#else
        std::mutex mutex;
        std::unique_lock<std::mutex> lk;
        std::condition_variable cond;
        std::vector<std::thread> threads;
#endif
        std::vector<checker<P>> _checkers;

    public:
        control_type(
                const typename checker<P>::method meth,
                const unsigned threads_num
        ) noexcept(false) : cond(), mutex() {
            _checkers = std::vector<checker<P>>(threads_num, checker<P>(this, meth));
#ifdef PTHREAD
            threads = std::vector<pthread_t>(threads_num);

            pthread_mutexattr_t attr;
            pthread_attr_t thread_attr;

            if (pthread_mutexattr_init(&attr) ||
                pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK) ||
                pthread_mutex_init(&mutex, &attr) ||
                pthread_cond_init(&cond, nullptr) ||
                pthread_attr_init(&thread_attr) ||
                pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED)) {
                throw std::runtime_error("pthread init failed");
            }

            pthread_mutex_lock(&mutex);

            for (unsigned i = 0; i < threads_num; ++i) {
                pthread_create(&threads[i], &thread_attr, &checker<P>::check, &_checkers[i]);
            }

            if (pthread_mutexattr_destroy(&attr) ||
                pthread_attr_destroy(&thread_attr)) {
                throw std::runtime_error("pthread destroy failed");
            }
#else
            lk = std::unique_lock<std::mutex>(mutex);
            threads = std::vector<std::thread>(threads_num);

            for (unsigned i = 0; i < threads_num; ++i) {
                threads[i] = std::thread(&checker<P>::check, &_checkers[i]);
                threads[i].detach();
            }
#endif
        }

#ifdef PTHREAD
        ~control_type() noexcept(false) {
            pthread_mutex_unlock(&mutex);
            if (pthread_cond_destroy(&cond) ||
                pthread_mutex_destroy(&mutex)) {
                throw std::runtime_error("pthread destroy failed");
            }
        }
#endif

        void lock() noexcept {
#ifdef PTHREAD
            pthread_mutex_lock(&mutex);
#else
            mutex.lock();
#endif
        }

        void unlock() noexcept {
#ifdef PTHREAD
            pthread_mutex_unlock(&mutex);
#else
            mutex.unlock();
#endif
        }

        void wait() noexcept {
#ifdef PTHREAD
            pthread_cond_wait(&cond, &mutex);
#else
            cond.wait(lk);
#endif
        }

        void signal() noexcept {
#ifdef PTHREAD
            pthread_cond_signal(&cond);
#else
            cond.notify_one();
#endif
        }

        std::vector<checker<P>> &checkers() noexcept {
            return _checkers;
        }
    };

private:
    polynomialgf<P> poly;

    control_type *ctrl;

    bool _busy;
    bool _terminate;
    result_type res;
    method _method;

public:
    explicit
    checker(control_type *, method = berlekamp) noexcept;

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

template<uint32_t P>
checker<P>::checker(control_type *ctrl, const method meth) noexcept :
        poly(), ctrl(ctrl), _busy(false), res({false, false}), _method(meth), _terminate(false) {}

template<uint32_t P>
polynomialgf<P> checker<P>::get() const {
    return poly;
}

template<uint32_t P>
void checker<P>::set(polynomialgf<P> val) {
    poly = val;
    _busy = true;
}

template<uint32_t P>
#ifdef PTHREAD
void *checker<P>::check(void *arg) noexcept {
    auto *c = static_cast<checker *>(arg);
#else
void checker<P>::check() noexcept {
    auto *c = this;
#endif
    while (true) {
        if (c->_terminate) { break; }
        if (!c->_busy) { continue; }
        switch (c->_method) {
            case method::rabin:
                c->res.irreducible = is_irreducible_rabin(c->poly);
                break;
            default: // method::berlekamp
                c->res.irreducible = is_irreducible_berlekamp(c->poly);
                break;
        }
        c->res.primitive = c->res.irreducible ? is_primitive(c->poly) : false;

        c->ctrl->lock();
        c->_busy = false;
        c->ctrl->signal();
        c->ctrl->unlock();
    }
#ifdef PTHREAD
    pthread_exit(nullptr);
#endif
}

template<uint32_t P>
[[nodiscard]]
bool checker<P>::busy() const noexcept {
    return _busy;
}

template<uint32_t P>
void checker<P>::terminate() noexcept {
    _terminate = true;
}

template<uint32_t P>
[[nodiscard]]
const typename checker<P>::result_type &checker<P>::result() const noexcept {
    return res;
}

#endif //IRRPOLY_CHECKER_HPP
