#ifndef IRRPOLY_CHECKER_HPP
#define IRRPOLY_CHECKER_HPP

#include <pthread.h>
#include "polynomialgf.hpp"

template<uint32_t P>
class checker {
    polynomialgf<P> poly;

    pthread_mutex_t *mutex;
    pthread_cond_t *cond;

    bool _busy;
public:
    struct result_type {
        bool irreducible;
        bool primitive;
    };

private:
    result_type res;

public:
    enum method {
        berlekamp = 0,
        rabin
    };

private:
    const method _method;

public:
    explicit
    checker(pthread_mutex_t *, pthread_cond_t *, method = berlekamp) noexcept;

    polynomialgf<P> get() const;

    void set(polynomialgf<P>);

    static
    void *check(void *arg) noexcept;

    [[nodiscard]]
    bool busy() const noexcept;

    [[nodiscard]]
    const result_type &result() const noexcept;
};

template<uint32_t P>
checker<P>::checker(pthread_mutex_t *mutex, pthread_cond_t *cond, const method meth) noexcept :
        poly(), mutex(mutex), cond(cond), _busy(false), res({false, false}), _method(meth) {}

template<uint32_t P>
polynomialgf<P> checker<P>::get() const {
    return poly;
}

template<uint32_t P>
void checker<P>::set(polynomialgf<P> val) {
    _busy = true;
    res = { false, false };
    poly = val;
}

template<uint32_t P>
void *checker<P>::check(void *arg) noexcept {
    auto *c = static_cast<checker *>(arg);

    switch (c->_method) {
        case method::rabin:
            c->res.irreducible = is_irreducible_rabin(c->poly);
            break;
        default: // method::berlekamp
            c->res.irreducible = is_irreducible_berlekamp(c->poly);
            break;
    }
    if (c->res.irreducible) {
        c->res.primitive = is_primitive(c->poly);
    }

    pthread_mutex_lock(c->mutex);
    c->_busy = false;
    pthread_cond_signal(c->cond);
    pthread_mutex_unlock(c->mutex);

    pthread_exit(nullptr);
}

template<uint32_t P>
bool checker<P>::busy() const noexcept {
    pthread_mutex_lock(mutex);
    auto b = _busy;
    pthread_mutex_unlock(mutex);
    return b;
}

template<uint32_t P>
const typename checker<P>::result_type &checker<P>::result() const noexcept {
    return res;
}

#endif //IRRPOLY_CHECKER_HPP
