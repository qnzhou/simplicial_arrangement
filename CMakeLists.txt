cmake_minimum_required(VERSION 3.14)

project(simplicial_arrangement LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(absl)
include(implicit_predicates)
include(spdlog)
include(nlohmann-json)
include(sanitizers)
include(cli11)
include(mshio)

option(SIMPLICIAL_ARRANGEMENT_UNIT_TESTS
    "Build unit tests for simplicial arrangement" OFF)
option(SIMPLICIAL_ARRANGEMENT_LOOKUP_TABLE_GEN
    "Lookup table generation mode (for developers only)" OFF)
option(SIMPLICIAL_ARRANGEMENT_DETERMINISTIC
    "Reqruire deterministic result (slightly slower)" OFF)
option(SIMPLICIAL_ARRANGEMENT_NON_ROBUST
    "Use non-robust predicates (for developers only)" OFF)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(simplicial_arrangement STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(simplicial_arrangement PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(simplicial_arrangement
    PUBLIC absl::numeric absl::flat_hash_map nlohmann_json::nlohmann_json
    PRIVATE implicit_predicates::implicit_predicates spdlog::spdlog)
target_compile_features(simplicial_arrangement PRIVATE cxx_std_17)
set(LOOKUP_TABLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lookup_table")
target_compile_definitions(simplicial_arrangement PUBLIC -DLOOKUP_TABLE_PATH="${LOOKUP_TABLE_PATH}")

if (SIMPLICIAL_ARRANGEMENT_LOOKUP_TABLE_GEN)
    target_compile_definitions(simplicial_arrangement PUBLIC -DLOOKUP_TABLE_GEN)
endif()
if (SIMPLICIAL_ARRANGEMENT_DETERMINISTIC)
    target_compile_definitions(simplicial_arrangement PRIVATE
        -DSIMPLICIAL_ARRANGEMENT_DETERMINISTIC)
endif()
if (SIMPLICIAL_ARRANGEMENT_NON_ROBUST)
    target_compile_definitions(simplicial_arrangement PRIVATE
        -DSIMPLICIAL_ARRANGEMENT_NON_ROBUST)
endif()
if (NOT MSVC)
    target_compile_options(simplicial_arrangement PRIVATE "-Wall" "-Werror")
endif()
add_library(simplicial_arrangement::simplicial_arrangement ALIAS
    simplicial_arrangement)

if (SIMPLICIAL_ARRANGEMENT_UNIT_TESTS)
    if (SIMPLICIAL_ARRANGEMENT_LOOKUP_TABLE_GEN)
        message(FATAL_ERROR "Lookup table generation mode cannot used for building unit tests.")
    endif()
    include(CTest)
    enable_testing()
    include(Catch2)

    file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
    add_executable(simplicial_arrangement_tests ${TEST_FILES})
    target_include_directories(simplicial_arrangement_tests PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/src/")
    target_link_libraries(simplicial_arrangement_tests PRIVATE
        simplicial_arrangement::simplicial_arrangement
        implicit_predicates::implicit_predicates Catch2::Catch2 spdlog::spdlog)
    target_compile_features(simplicial_arrangement_tests PRIVATE cxx_std_17)
    if (NOT MSVC)
        target_compile_options(simplicial_arrangement_tests PRIVATE "-Wall" "-Werror")
    else()
        target_compile_definitions(simplicial_arrangement_tests PRIVATE "-D_USE_MATH_DEFINES")
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
    if (SIMPLICIAL_ARRANGEMENT_NON_ROBUST)
        target_compile_definitions(simplicial_arrangement_tests PRIVATE
            -DSIMPLICIAL_ARRANGEMENT_NON_ROBUST)
    endif()
endif()

if (SIMPLICIAL_ARRANGEMENT_LOOKUP_TABLE_GEN)
    add_executable(generate_material_interface_lookup_table
        "app/generate_material_interface_lookup_table.cpp")
    target_link_libraries(generate_material_interface_lookup_table PRIVATE
        simplicial_arrangement::simplicial_arrangement nlohmann_json::nlohmann_json
        spdlog::spdlog)
    target_compile_features(generate_material_interface_lookup_table PRIVATE cxx_std_17)

    add_executable(generate_simplicial_arrangement_lookup_table
        "app/generate_simplicial_arrangement_lookup_table.cpp")
    target_link_libraries(generate_simplicial_arrangement_lookup_table PRIVATE
        simplicial_arrangement::simplicial_arrangement nlohmann_json::nlohmann_json
        spdlog::spdlog)
    target_compile_features(generate_simplicial_arrangement_lookup_table PRIVATE cxx_std_17)
endif()

add_executable(simplex_case_study "app/simplex_case_study.cpp")
target_link_libraries(simplex_case_study PRIVATE
    simplicial_arrangement::simplicial_arrangement mshio::mshio CLI11::CLI11
    nlohmann_json::nlohmann_json spdlog::spdlog)
target_compile_features(simplex_case_study PRIVATE cxx_std_17)
