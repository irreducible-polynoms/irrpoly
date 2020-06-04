cmake_minimum_required(VERSION 3.8)

include(currdir) # save current directory name in CURRDIR variable
project("${PROJECT_NAME}.${CURRDIR}")

include(sources)
