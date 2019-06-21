cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(MMseqs_benchmark)
set(HAVE_FFINDEX 1 CACHE BOOL "Have FFindex")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

include_directories(lib/kseq)

add_executable(evaluate_results EvaluateResults.cpp EvaluateResults.h)
add_executable(evaluate_pos_results EvaluateResultsPos.cpp EvaluateResults.h)
add_executable(test_main TestMain.cpp)

if (HAVE_FFINDEX)
    add_definitions(-DHAVE_FFINDEX)
    add_executable(create_dataset CreateDataset.cpp)
    add_executable(create_dataset_pos CreateDatasetPos.cpp)

    add_subdirectory(lib/ffindex)
    include_directories(lib/ffindex/src)

    target_link_libraries(create_dataset ffindex)
    target_link_libraries(create_dataset_pos ffindex)
    target_link_libraries(evaluate_results ffindex)
    target_link_libraries(evaluate_pos_results ffindex)
endif()


install(TARGETS 
    evaluate_results
    evaluate_pos_results
    test_main
    DESTINATION bin)
