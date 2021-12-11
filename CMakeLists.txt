cmake_minimum_required(VERSION 3.14)

project(simplicial_arrangement LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(absl)
include(implicit_predicates)
include(spdlog)
include(nlohmann-json)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(simplicial_arrangement STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(simplicial_arrangement PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(simplicial_arrangement
    PUBLIC absl::numeric absl::flat_hash_map nlohmann_json::nlohmann_json
    PRIVATE implicit_predicates::implicit_predicates spdlog::spdlog)
target_compile_features(simplicial_arrangement PRIVATE cxx_std_17)
if (NOT MSVC)
    target_compile_options(simplicial_arrangement PRIVATE "-Wall" "-Werror")
endif()
add_library(simplicial_arrangement::simplicial_arrangement ALIAS
    simplicial_arrangement)

option(SIMPLICIAL_ARRANGEMENT_UNIT_TESTS "Build unit tests for simplicial arrangement" OFF)
if (SIMPLICIAL_ARRANGEMENT_UNIT_TESTS)
    include(CTest)
    enable_testing()
    include(Catch2)

    file(GLOB TEST_FILES "${CMAKE_SOURCE_DIR}/tests/*.cpp")
    add_executable(simplicial_arrangement_tests ${TEST_FILES})
    target_link_libraries(simplicial_arrangement_tests PRIVATE
        simplicial_arrangement::simplicial_arrangement
        implicit_predicates::implicit_predicates Catch2::Catch2)
    target_compile_features(simplicial_arrangement_tests PRIVATE cxx_std_17)
    if (NOT MSVC)
        target_compile_options(simplicial_arrangement_tests PRIVATE "-Wall" "-Werror")
    endif()
    target_compile_definitions(simplicial_arrangement_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
    catch_discover_tests(simplicial_arrangement_tests)

    option(SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE "Track unit test coverage" OFF)
    if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
        if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
            target_compile_options(simplicial_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
            target_link_options(simplicial_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
            target_compile_options(simplicial_arrangement_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
            target_link_options(simplicial_arrangement_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        endif()
    endif()
endif()


set(LOOKUP_TABLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lookup_table")
target_compile_definitions(simplicial_arrangement PUBLIC -DLOOKUP_TABLE_PATH="${LOOKUP_TABLE_PATH}")

add_executable(test_lookup_table "app/test_lookup_table.cpp")
target_link_libraries(test_lookup_table simplicial_arrangement::simplicial_arrangement implicit_predicates::implicit_predicates nlohmann_json::nlohmann_json)

add_executable(impl_arrangement "app/implicit_arrangement.cpp")
target_link_libraries(impl_arrangement simplicial_arrangement::simplicial_arrangement implicit_predicates::implicit_predicates nlohmann_json::nlohmann_json)