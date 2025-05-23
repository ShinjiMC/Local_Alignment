cmake_minimum_required(VERSION 3.10)
project(AlignmentGlobal)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_executable(align
    main.cpp
    src/localAligner.cpp
    src/alignmentIO.cpp
)


# ===============================
# Google Test: configuration test
# ===============================
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set(TEST_SOURCES
    test.cpp
    src/localAligner.cpp
    src/alignmentIO.cpp
)

add_executable(main_test ${TEST_SOURCES})

target_link_libraries(main_test
    ${GTEST_LIBRARIES}
    pthread
)

target_link_options(main_test PRIVATE "-Wl,-rpath,/usr/lib64")

enable_testing()
add_test(NAME main_test COMMAND main_test)