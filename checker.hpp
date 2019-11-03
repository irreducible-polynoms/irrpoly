#ifndef IRRPOLY_CHECKER_HPP
#define IRRPOLY_CHECKER_HPP

#include <pthread.h>
#include "polynomialgf.hpp"

template<uint32_t P>
class checker {
public:
    struct control_type {
        pthread_mutex_t mutex{};
        pthread_cond_t cond{};
        std::vector<checker<P>> checkers;
        std::vector<pthread_t> threads;
    };

    struct result_type {
        bool irreducible;
        bool primitive;
    };

    enum method {
        berlekamp = 0,
        rabin
    };

private:
    polynomialgf<P> poly;

    pthread_mutex_t *mutex;
    pthread_cond_t *cond;

    bool _busy;
    bool _terminate;
    result_type res;
    method _method;

public:
    explicit
    checker(pthread_mutex_t *, pthread_cond_t *, method = berlekamp) noexcept;

    polynomialgf<P> get() const;

    void set(polynomialgf<P>);

    [[nodiscard]]
    bool busy() const noexcept;

    void terminate() noexcept;

    [[nodiscard]]
    const result_type &result() const noexcept;

    static
    void *check(void *arg) noexcept;

    static
    int init(control_type &, method, unsigned) noexcept;

    static
    int destroy(control_type &) noexcept;
};

template<uint32_t P>
checker<P>::checker(pthread_mutex_t *mutex, pthread_cond_t *cond, const method meth) noexcept :
        poly(), mutex(mutex), cond(cond), _busy(false), res({false, false}), _method(meth), _terminate(false) {}

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
void *checker<P>::check(void *arg) noexcept {
    auto *c = static_cast<checker *>(arg);

    while (true) {
        if (c->_terminate) { pthread_exit(nullptr); }
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

        pthread_mutex_lock(c->mutex);
        c->_busy = false;
        pthread_cond_signal(c->cond);
        pthread_mutex_unlock(c->mutex);
    }
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

template<uint32_t P>
int checker<P>::init(
        typename checker<P>::control_type &ctrl,
        const typename checker<P>::method meth,
        const unsigned threads_num
) noexcept {
    ctrl.checkers = std::vector<checker<P>>(threads_num, checker<P>(&ctrl.mutex, &ctrl.cond, meth));
    ctrl.threads = std::vector<pthread_t>(threads_num);

    pthread_mutexattr_t attr;
    pthread_attr_t thread_attr;

    if (pthread_mutexattr_init(&attr) ||
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK) ||
        pthread_mutex_init(&ctrl.mutex, &attr) ||
        pthread_cond_init(&ctrl.cond, nullptr) ||
        pthread_attr_init(&thread_attr) ||
        pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED)) {
        return -1;
    }

    for (unsigned i = 0; i < threads_num; ++i) {
        pthread_create(&ctrl.threads[i], &thread_attr, &checker<P>::check, &ctrl.checkers[i]);
    }

    if (pthread_mutexattr_destroy(&attr) ||
        pthread_attr_destroy(&thread_attr)) {
        return -1;
    }
    return 0;
}

template<uint32_t P>
int checker<P>::destroy(checker::control_type &ctrl) noexcept {
    if (pthread_cond_destroy(&ctrl.cond) ||
        pthread_mutex_destroy(&ctrl.mutex)) {
        return -1;
    }
    return 0;
}

#endif //IRRPOLY_CHECKER_HPP
