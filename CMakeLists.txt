cmake_minimum_required(VERSION 2.6)
#set(CMAKE_CXX_COMPILER "/usr/bin/g++")	# to avoid using clang on macbook
#set(CMAKE_C_COMPILER "/usr/bin/gcc")	# ,,

set(CMAKE_INSTALL_PREFIX "$ENV{HOME}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

project(mclattice2)


include_directories("./src/")


## Standard
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -DNDEBUG")
## Debug
#set(CMAKE_CXX_FLAGS "-std=c++11 -O0 -DDEBUG -g")

#set(CMAKE_CXX_FLAGS "-std=c++11 -O0 -DDEBUG -g")
#set(CMAKE_CXX_FLAGS "-std=c++11")
# set (CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -DDEBUG" ${CMAKE_CXX_FLAGS_DEBUG})

find_package(HDF5 REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

include_directories(${HDF5_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
set(_hdf5_libs hdf5 hdf5_cpp)


file(GLOB mclattice2_SRC
  src/*.cpp
  src/*/*.cpp
)
add_executable(mclattice2 ${mclattice2_SRC} ${sources})


add_library(enhance STATIC lib/enhance.cpp)


target_link_libraries(mclattice2 enhance)
target_link_libraries(mclattice2 ${_hdf5_libs})
target_link_libraries(mclattice2 ${Boost_LIBRARIES})

install(TARGETS mclattice2 RUNTIME DESTINATION bin)

