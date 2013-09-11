cmake_minimum_required( VERSION 2.6 )
set(CMAKE_CXX_FLAGS "-Wall -O2")

project( genetic-ifs )
enable_testing ()

include_directories( include )
#include_directories( ext/gtest/include )
add_subdirectory( src )
#add_subdirectory( ext/gtest )
#add_subdirectory( test )
