cmake_minimum_required(VERSION 4.0.1)
set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_BUILD_TYPE Debug)

project(TriangleIntersection VERSION 1.0.0)

add_library(TriangleIntersection)

target_compile_features(TriangleIntersection PUBLIC cxx_std_11)

target_sources(TriangleIntersection
    PRIVATE
        src/intersection.cpp
        src/misc.cpp
        src/intersection_private.hpp

    PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
            inc
        FILES
            inc/intersection.hpp
)


add_executable(TriangleIntersectionTest)
target_sources(TriangleIntersectionTest
    PRIVATE
    test/test.cpp
)
target_compile_features(TriangleIntersection PUBLIC cxx_std_17)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
)

FetchContent_MakeAvailable(googletest)

enable_testing()

target_link_libraries(TriangleIntersectionTest GTest::gtest_main TriangleIntersection)

include(GoogleTest)
gtest_discover_tests(TriangleIntersectionTest)
