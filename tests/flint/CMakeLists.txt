cmake_minimum_required(VERSION 3.8)

include(currdir) # save current directory name in CURRDIR variable
project("${PROJECT_NAME}.${CURRDIR}")

# GMP
include_directories(/usr/local/Cellar/gmp/6.2.0/include)
link_directories(/usr/local/Cellar/gmp/6.2.0/lib)

# MPFR
include_directories(/usr/local/Cellar/mpfr/4.0.2/include)
link_directories(/usr/local/Cellar/mpfr/4.0.2/lib)

# FLINT
include_directories(/usr/local/Cellar/flint/2.6.0/include)
link_directories(/usr/local/Cellar/flint/2.6.0/lib)

#include(sources)
#target_link_libraries("${PROJECT_NAME}" gmp mpfr flint)
