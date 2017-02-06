cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(MMseqs_benchmark)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES CreateDataset.cpp)
add_executable(create_dataset ${SOURCE_FILES})

# libs
add_subdirectory(lib/ffindex)
include_directories(${PROJECT_SOURCE_DIR}/lib/ffindex/src)
include_directories(${PROJECT_SOURCE_DIR}/lib/kseq/)

target_link_libraries(create_dataset ffindex)

add_executable(evaluate_results EvaluateResults.cpp EvaluateResults.h)
target_link_libraries(evaluate_results ffindex)

add_executable(test_main TestMain.cpp)
