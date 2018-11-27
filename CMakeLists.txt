#cmake_minimum_required(VERSION 3.6)
project(fastcky)

if($ENV{CXX})
  set(CMAKE_CXX_COMPILER $ENV{CXX} CACHE FILEPATH "CXX Compiler")
endif()
if($ENV{CC})
  set(CMAKE_CC_COMPILER $ENV{CC} CACHE FILEPATH "CC Compiler")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lprofiler")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltcmalloc -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")

set(BEAMCKYPAR_SOURCE
        beamckypar.cpp
        beamckypar.h
        fastcky_w.h
        utility.h
        cxxopts.hpp)

add_executable(linearfold ${BEAMCKYPAR_SOURCE})
