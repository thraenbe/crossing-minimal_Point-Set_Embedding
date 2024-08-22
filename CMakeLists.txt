cmake_minimum_required(VERSION 3.13)
project(PointSetCrossMin)

set(CMAKE_CXX_STANDARD 20)

set(BOOST_ROOT "~/boost-ver/boost_1_85_0")
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(LIBXLSXWRITER_INCLUDE_DIRS "/usr/local/include")
set(LIBXLSXWRITER_LIBRARIES "/usr/lib/x86_64-linux-gnu/libxlsxwriter.so")

# Include the directories for the library
include_directories(${LIBXLSXWRITER_INCLUDE_DIRS})



add_executable(PointSetCrossMin
        main.cpp
        JSONParser.cpp
	    JSONParser.hpp
        Geometry.hpp
        json.hpp
        Graph.hpp
        matching.hpp
        Node.hpp
        )

target_include_directories(PointSetCrossMin PRIVATE ogdf/include ${Boost_INCLUDE_DIRS} ${LIBXLSXWRITER_LIBRARIES})
find_library(OGDF OGDF HINTS ogdf ./ogdf/Release)
find_library(COIN COIN HINTS ogdf ./ogdf/Release)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
    target_link_libraries(PointSetCrossMin PRIVATE ${OGDF} ${COIN} stdc++fs ${Boost_LIBRARIES}  ${LIBXLSXWRITER_LIBRARIES})
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    set(CMAKE_CXX_FLAGS "/Ox /EHsc /openmp")
    target_link_libraries(PointSetCrossMin PRIVATE ${OGDF} ${COIN} ${Boost_LIBRARIES}  ${LIBXLSXWRITER_LIBRARIES})
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
