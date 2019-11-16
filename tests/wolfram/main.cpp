#include <iostream>
#include <thread>

#include "irrpoly/polynomialgf.hpp"

using namespace irrpoly;

/// Функция, формирующая данные для тестов. Результаты скопировать в файл test_input.
void test() {
    {
        const uint32_t P = 2;
        for (size_t i = 4; i < 24; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random < P > (i);
                std::cout << (is_irreducible_rabin(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 3;
        for (size_t i = 4; i < 18; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random < P > (i);
                std::cout << (is_irreducible_rabin(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 5;
        for (size_t i = 4; i < 14; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random < P > (i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 7;
        for (size_t i = 4; i < 10; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random < P > (i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 11;
        for (size_t i = 4; i < 6; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random < P > (i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
    {
        const uint32_t P = 13;
        for (size_t i = 4; i < 4; ++i) {
            for (size_t j = 0; j < i / 2; ++j) {
                polynomialgf<P> p = random < P > (i);
                std::cout << (is_irreducible_berlekamp(p) ? '1' : '0') << " "
                          << (is_primitive_definition(p) ? '1' : '0') << " " << P << " "
                          << p << std::endl;
            }
        }
    }
}

int main() {
    std::cout << "TEST DATA" << std::endl;
    test();
    return 0;
}
