# Irreducibility and primitivity tests for polynomials over finite fields

## Brief
This project implements several irreducibility and primitivity tests
for polynomials in finite fields. Functions in this library are currently
faster then similar function in Wolfram and Mathlab, but slower than
functions from PARI/GP and FLINT. This library has a guarantee of
NO UNDEFINED BEHAVIOUR by design.

Currently implemented tests:
1. Berlekamp's irreducibility test
2. Rabin's irreducibility test
3. Ben-Or's irreducibility test
4. Primitivity test

Recommended functions to use:
- `is_irreducible` for irreducibility test
- `is_primitive` for primitivity test

## Requirements
- C++ 17 support
- Cmake 3.8 or higher

## Limitations
- Only PRIME fields are supported
- Multithreading is currently slow and not recommended for use

## Usage
1. Download the newest sources from
    [Releases](https://github.com/irreducible-polynoms/irrpoly/releases)
    and unpack them to you root project folder (so that path to sources would be `./irrpoly`).
    Another way is using Git Submodules with command
    ```git
    git submodule add https://github.com/irreducible-polynoms/irrpoly
   ```
2. Check that your `CMAKE_CXX_STANDARD` is 17 or greater.
3. Add the following to your `CMakeLists.txt`:
    ```C++
    add_library(irrpoly INTERFACE)
    target_include_directories(irrpoly INTERFACE irrpoly/include)
    add_dependencies(irrpoly Threads::Threads)
    ```
    If you are not planning to use miltithread versions of algorithms be free
    to remove the last line. Now to link with `irrpoly` library you need to add
    ```C++
    target_link_libraries("${PROJECT_NAME}" irrpoly)
    ```
    after `add_executable` (replace `${PROJECT_NAME}` with your target's name if it differs).
4. Add `#include <irrpoly.h>` to use the library. Remember that everything is dived into
    `namespace irrpoly`.
5. If you want to enable field consistency checks in Release configuration - place
    `#define IRRPOLY_RELEASE_CHECKED` before `#include <irrpoly.h>`. By default this checks
    are enabled only for Debug configuration to speed up Release.

## Contents
- `gf` – represents Galois field, to create new instance of `gf` use `make_gf` function
- `gfn` – represents a number in Galois field
- `gfpoly` – represents a polynomial with coefficients from Galois field
- `gfcheck` – contains checks implementations and some helpers (`gcd`, `derivative`)
- `pipeline` – contains multithread pipeline for performing polynomial checks
    (examples of usage provided [here](examples)), everything dived in namespace
    `irrpoly::multithread`
- `nn` – redistributed `dropbox::oxygen::nn` class
    ([original source](https://github.com/dropbox/nn))

## Documentation
Auto-generated documentation is placed [here](docs) and could be accessed as
[GitHub Pages](https://irreducible-polynoms.github.io/irrpoly/html/).

## For developers
- Building `tests` and `examples` executables if possible using [Makefile](Makefile)
    command `make build` of from any IDE supporting Cmake (JetBrains CLion, Visual Studio, etc.).
    Binaries are placed near sources in folders `bin/debug` and `bin/release`.
- Unit-tests and benchmarks are written using [Catch2](https://github.com/catchorg/Catch2)
- If PVS-Studio is installed - targets `*.analyze` are automatically added for all
    examples, use them for performing static analysis for possible vulnerabilities.
    When adding new examples preserve the following comment at the top of all `*.cpp` files:
    ```C++
    // This is an open source non-commercial project. Dear PVS-Studio, please check it.
    // PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
    ```
- Updating documentation requires Doxygen and could be performed with `make docs`.
- When bumping version number remember to change it both in [CMakeLists.txt](CMakeLists.txt) 
    and in [Doxyfile](docs/Doxyfile), and make a proper tag (or GitHub Release).
- If you are going to change check methods - verify their correctness using method
    described [here](tests/wolfram/README.md).

## License
The `irrpoly` library is distributed under MIT license,
some parts of code are from `boost` library originally distributed under Boost license
(redistributed under MIT license),
`nn` class is originally distributed under Apache license (redistributed under MIT license).
`Catch2` is not redistributed as it is only used for tests.

## TODO
- Create a script for generating single-header-version of library and 
    add it to [Makefile](Makefile)'s `make build` command (check is it allowed
    to redistribute `nn` in such form).
- Replace slow Berlekamp's matrix rank calculation algorithm in
    `is_irreducible_berlekamp` method with faster one.
- Rewrite `pipeline` to make it faster (keep in mind that for most polynomials
    check functions will return `false` very quickly). C++ 20 coroutines could
    be used here.
- Add support of COMPOSITE fields.
- Implement Distinct Degree Factorisation algorithm for irreducibility check.
- Research the possibility to use Discrete Fourier Transform to speed up
    polynomial multiplication and division methods.
- Check [FLINT](http://www.flintlib.org/) sources and find out the way Rabin's
    irreducibility test is implemented there (there is some analogue of `x_pow_mod`
    function from this library), also inspect other methods as they are
    significantly faster then (`nmod_poly_is_irreducible` stands for Berlekamp-Zassenhaus
    algorithm, `nmod_poly_is_irreducible_ddf` stands for Distinct Degree Factorisation
    algorithm, `nmod_poly_is_irreducible_rabin` stands for Rabin's irreducibility test)

## Literature
1. [Wikipedia page](https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields)
    on factorization of polynomials over finite fields.
2. Davison, Joosten, Thiemann, Yamada. A Formalization of Berlekamp’s Factorization Algorithm.
3. Panario, Pittel, Richmond, Viola. Analysis of Rabin's irreducibility test for polynomials
    over finite Fields.
4. Gao, Panario. Tests and constructions of irreducible polynomials over finite fields.
5. Hansen, Mullen. Primitive polynomials over finite fields.
