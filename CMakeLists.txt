cmake_minimum_required(VERSION 3.13)
project(PointSetCrossMin)

set(CMAKE_CXX_STANDARD 20)

add_executable(PointSetCrossMin
        main.cpp
        JSONParser.cpp
	    JSONParser.hpp
        Geometry.hpp
        json.hpp
        Graph.hpp
        )

target_include_directories(PointSetCrossMin PRIVATE ogdf/include)
find_library(OGDF OGDF HINTS ogdf ./ogdf/Release)
find_library(COIN COIN HINTS ogdf ./ogdf/Release)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
    target_link_libraries(PointSetCrossMin PRIVATE ${OGDF} ${COIN} stdc++fs)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    set(CMAKE_CXX_FLAGS "/Ox /EHsc /openmp")
    target_link_libraries(PointSetCrossMin PRIVATE ${OGDF} ${COIN})
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
