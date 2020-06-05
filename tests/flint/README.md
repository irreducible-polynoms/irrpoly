# Check methods' correctness verification using FLINT
This method requires FLINT library and it's dependencies GMP (MPIR)
and MPFR to be installed. It you've installed them - update their
`include` and `lib` paths in [CMakeLists.txt](CMakeLists.txt)
and uncomment lines:
```CMake
#include(sources)
#target_link_libraries("${PROJECT_NAME}" gmp mpfr flint)
```
If you haven't installed required libraries and uncommented this
lines - build errors will occur and whole CMake project won't be
build.

Known issues: FLINT returns `true` for constants, while
this project and Wolfram returns `false`.